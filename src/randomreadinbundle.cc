/**   LICENCE
* Copyright (c) 2020 Genome Research Ltd.
* 
* Author: Cancer Genome Project <cgphelp@sanger.ac.uk>
* 
* This file is part of NanoSeq.
* 
* NanoSeq is free software: you can redistribute it and/or modify it under
* the terms of the GNU Affero General Public License as published by the Free
* Software Foundation; either version 3 of the License, or (at your option) any
* later version.
* 
* This program is distributed in the hope that it will be useful, but WITHOUT
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
* FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
* details.
* 
* You should have received a copy of the GNU Affero General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

/*

Run:
  ./randomreadinbundle -I in.bam -O out.bam -m 1

*/

#include "randomreadinbundle.h"

std::vector<std::string> tokenize(std::string str, char delimiter) {
  std::istringstream iss(str);
  std::vector<std::string> tokens;
  std::string token;
  while (std::getline(iss, token, delimiter)) {
    tokens.push_back(token);
  }
  return tokens;
}

void RandomlySelectOneReadPerBundle::LoadFiles() {
  this->in = hts_open(this->infile, "r");
  if (this->in == 0) {
    std::stringstream er;
    er << "Fail to open input BAM/CRAM file ";
    er << this->infile;
    er << std::endl;
    throw std::invalid_argument(er.str());
  }
  this->head = sam_hdr_read(this->in);
  if (this->head == NULL || this->head->n_targets == 0) {
    std::stringstream er;
    er << "Error: BAM/CRAM file does not have header.";
    er << std::endl;
    throw std::runtime_error(er.str());
  }
  if (BamIsCorrectlyPreprocessed(head) == false) {
      std::stringstream er;
      er << "Error : bam is not properly preprocessed. (missing RG: tag from bamreadbundles)";
      er << std::endl;
      throw std::runtime_error(er.str());
    }
  this->out = hts_open(this->outfile, "wb");
  if (!this->out) {
    std::stringstream er;
    er << "Error: failed to open ";
    er << this->outfile;
    er << " for output.";
    er << std::endl;
    throw std::invalid_argument(er.str());
  }
}

bool RandomlySelectOneReadPerBundle::BamIsCorrectlyPreprocessed(bam_hdr_t *head ) {
  std::vector<std::string> tokens = tokenize(head->text, '\n');
  for (int j = 0; j < tokens.size(); j++) {
    if (tokens[j].rfind("@PG", 0) == 0) {
      std::vector<std::string> subtokens = tokenize(tokens[j], '\t');
      for (int k = 0; k < subtokens.size(); k++) {
        if (subtokens[k].rfind("ID", 0) == 0) {
          if (subtokens[k].rfind("ID:bamaddreadbundles", 0) == 0) {
            return true;
          }
        }
      }
    }
  }
  return false;
}

void RandomlySelectOneReadPerBundle::UpdateHeader(int argc, char **argv) {
  // command line argument
  std::stringstream ss;
  int i;
  for (i = 0; i < (argc-1); ++i) {
    ss << argv[i];
    ss << " ";
  }
  ss << argv[i];
  sam_hdr_t *sh = sam_hdr_parse(this->head->l_text, this->head->text);
  char pg[]   = "randomreadinbundle";
  sam_hdr_add_pg(sh, pg, "VN", "1.0", "CL", ss.str().c_str(), NULL);
  free(this->head->text);
  sam_hdr_rebuild(sh);
  this->head->text   = strdup(sam_hdr_str(sh));
  this->head->l_text = sam_hdr_length(sh);
  sam_hdr_destroy(sh);
  if (sam_hdr_write(this->out, this->head) < 0) {
      std::stringstream er;
      er << "Failed to write header";
      er << std::endl;
      throw std::runtime_error(er.str());
  }
}

int RandomlySelectOneReadPerBundle::ReadStrand(bam1_t* b) {
  if (b->core.flag & BAM_FMREVERSE) {
    return 0;
  } else if (b->core.flag & BAM_FREVERSE) {
    return 1;
  } else {
    throw std::invalid_argument("Invalid strand");
  }
}


void RandomlySelectOneReadPerBundle::WriteOut(bam1_t* b) {
  if ((sam_write1(this->out, this->head, b) < 0)) {
    std::stringstream er;
    er << "Error: failed to write record.";
    er << std::endl;
    throw std::runtime_error(er.str());
  }
}


void RandomlySelectOneReadPerBundle::SelectReads() {
  bam1_t *b       = bam_init1();
  std::map<std::string, std::map<int, int>> seen;
  int last_chr = 0;
  int i_chr = 0; 
  while (sam_read1(this->in, this->head, b) >= 0) {
    i_chr = b->core.tid;
    if ( i_chr != last_chr ) {
      seen.clear();
      last_chr = i_chr;
    }
    std::string rb = std::string(bam_aux2Z(bam_aux_get(b, "RB")));
    int strand     = RandomlySelectOneReadPerBundle::ReadStrand(b);
    if (seen.count(rb) == 0) {
      seen[rb][0] = 0;
      seen[rb][1] = 0;
      seen[rb][strand]++;
    } else {
      seen[rb][strand]++;
    }

    if (seen[rb][strand] == this->min_reads) {
	    b->core.flag &= ~BAM_FDUP;
      RandomlySelectOneReadPerBundle::WriteOut(b);
    }
  }
  seen.clear();
  bam_destroy1(b);

  bam_hdr_destroy(this->head);
  hts_close(this->in);
  if (sam_close(this->out) < 0) {
    std::stringstream er;
    er << "Error closing output file.";
    er << std::endl;
    throw std::invalid_argument(er.str());
  }
}


void Usage() {
  fprintf(stderr, "\nUsage:\n");
  fprintf(stderr, "\t-I\tInput BAM/CRAM file name\n");
  fprintf(stderr, "\t-O\tOutput BAM file name\n");
  fprintf(stderr, "\t-m\tMinimum number of reads in read-bundle\n");
  fprintf(stderr, "\t-h\tHelp\n");
}


int main(int argc, char **argv) {
  RandomlySelectOneReadPerBundle rs;
  rs.infile = NULL;
  rs.outfile = NULL;
  rs.min_reads = 1;
  int opt = 0;
  while ((opt = getopt(argc, argv, "I:O:m:h")) >= 0) {
    switch (opt) {
      case 'I':
        rs.infile = optarg;
        break;
      case 'O':
        rs.outfile = optarg;
        break;
      case 'm':
        rs.min_reads = std::stoi(optarg);
        break;
      case 'h':
        Usage();
        exit(0);
      default:
        break;
    }
  }
  if (rs.infile == NULL) {
    std::stringstream er;
    er << "Error: no input file specified";
    throw std::invalid_argument(er.str());
  }
  if (rs.outfile == NULL) {
    std::stringstream er;
    er << "Error: no output file specified";
    throw std::invalid_argument(er.str());
  }
  rs.LoadFiles();
  rs.UpdateHeader(argc, argv);
  rs.SelectReads();
  return 0;
}

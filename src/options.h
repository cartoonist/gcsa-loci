/**
 *    @file  options.h
 *   @brief  Options class definition.
 *
 *  Data structure for storing command-line option values.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Thu Aug 03, 2017  04:44
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef OPTIONS_H__
#define OPTIONS_H__

#include <string>

namespace gloci {
  typedef struct {
    std::string vg_filename;
    std::string gcsa_filename;
    std::string output_filename;
    unsigned int seed_len;
  } Options;
}

#endif  // GLOCI_OPTIONS_H__

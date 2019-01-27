/**
 *    @file  release.h
 *   @brief  Package release information.
 *
 *  The meta data about the package and its release information.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Thu Aug 03, 2017  17:06
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef RELEASE_H__
#define RELEASE_H__

#ifndef GIT_VERSION
#define GIT_VERSION "NONE"
#endif

#ifndef UPDATE_DATE
#define UPDATE_DATE "NONE"
#endif

namespace gloci{
  namespace release {
      const char* const version = GIT_VERSION;                 /**< @brief Version number. */
      const char* const name = "gcsa_loci";              /**< @brief Package name */
      const char* const short_desc = "GCSA2 loci coverage calculator";  /**< @brief Short description. */
      /** @brief Long description. */
      const char* const desc = "Calculate loci coverage of a GCSA2 index";
    }  /* -----  end of namespace package  ----- */
}

#endif  // RELEASE_H__

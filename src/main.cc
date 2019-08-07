/**
 *    @file  main.cc
 *   @brief  gcsa_loci main program.
 *
 *  Find uncovered loci in a GCSA2 index.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Thu Aug 03, 2017  04:37
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <cstdlib>
#include <csignal>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <seqan/arg_parse.h>
#include <sdsl/bit_vectors.hpp>
#include <gcsa/gcsa.h>
#include "grem/vargraph.h"
#include "grem/utils.h"

#include "timer.h"
#include "options.h"
#include "release.h"


using namespace gloci;

  seqan::ArgumentParser::ParseResult
parse_args( Options& options, int argc, char* argv[] );

  void
setup_argparser( seqan::ArgumentParser& parser );

  inline void
get_option_values( Options& options, seqan::ArgumentParser& parser );

  void
query_loci( std::string& vg_name, std::string& gcsa_name, unsigned int seed_len,
    std::string& output_name );

  void
signal_handler( int signal );

std::size_t done_idx = 0;
std::size_t total_no = 0;
std::size_t uncovered_loci = 0;


  int
main( int argc, char* argv[] )
{
  // Parse the command line.
  Options options;
  auto res = parse_args( options, argc, argv );
  // If parsing was not successful then exit with code 1 if there were errors.
  // Otherwise, exit with code 0 (e.g. help was printed).
  if (res != seqan::ArgumentParser::PARSE_OK)
    return res == seqan::ArgumentParser::PARSE_ERROR;

  /* Install signal handler */
  std::signal( SIGUSR1, signal_handler );

  query_loci( options.vg_filename, options.gcsa_filename, options.seed_len,
      options.output_filename );

  return EXIT_SUCCESS;
}


  void
signal_handler( int )
{
  std::cout << "Queried " << ::done_idx << " out of " << ::total_no
            << " with " << ::uncovered_loci << " uncovered loci in "
            << Timer<>::get_lap_str( "query" ) << ": "
            << ::done_idx * 100 / total_no << "% done." << std::endl;
}


  void
query_loci( std::string& vg_name, std::string& gcsa_name, unsigned int seed_len,
    std::string& /* output_name */ )
{
  typedef std::make_unsigned< grem::VarGraph::offset_type >::type uoffset_type;

  //std::ofstream output_file( output_name, std::ofstream::out );
  std::ifstream vg_file( vg_name, std::ifstream::in | std::ifstream::binary );
  if ( !vg_file ) {
    throw std::runtime_error("could not open file '" + vg_name + "'" );
  }
  std::ifstream gcsa_file( gcsa_name, std::ifstream::in | std::ifstream::binary );
  if ( !gcsa_file ) {
    throw std::runtime_error("could not open file '" + gcsa_name + "'" );
  }

  grem::VarGraph vargraph;
  gcsa::GCSA index;
  std::vector< std::string > sequences;
  std::vector< std::string > patterns;
  std::vector< gcsa::node_type > results;

  std::cout << "Loading graph..." << std::endl;
  if ( grem::ends_with( vg_name, ".vg" ) ) {
    vargraph.from_stream( vg_file );
  }
  else {
    vargraph.load( vg_file );
  }
  std::cout << "Loading GCSA index..." << std::endl;
  index.load( gcsa_file );

  total_no = vargraph.get_total_nof_loci();
  std::cout << "Querying loci..." << std::endl;
  {
    auto timer = Timer<>( "query" );
    seqan::Iterator< grem::VarGraph, grem::Backtracker >::Type bt_itr( vargraph );
    grem::Path< grem::VarGraph > trav_path( &vargraph );
    sdsl::bit_vector bv( vargraph.get_max_node_len(), 0 );
    for ( grem::VarGraph::rank_type rank = 1; rank <= vargraph.max_node_rank(); ++rank ) {
      grem::VarGraph::nodeid_type id = vargraph.rank_to_id( rank );
      auto label = vargraph.node_sequence( id );
      uoffset_type offset = ( label.size() < seed_len ) ? 0 : label.size() - seed_len + 1;

      gcsa::range_type range;
      // Query all k-mers inside the node.
      auto end = label.begin() + offset;
      for ( auto it = label.begin(); it != end; ++it ) {
        range = index.find( it, it + seed_len );
        // :TODO:Wed Aug 07 11:09:\@cartoonist:
        // Getting non-empty range doesn't necessarily mean that the locus is covered.
        // The locations corresponding to the range should also be verified.
        if ( gcsa::Range::empty( range ) ) ++uncovered_loci;
      }
      // Query all k-mers spanning between this node and other ones.
      go_begin( bt_itr, id );
      while ( !at_end( bt_itr ) && offset != label.size() ) {
        extend_to_k( trav_path, bt_itr, label.size() - 1 + seed_len );
        if ( trav_path.get_sequence_len() >= seed_len ) {
          auto trav_seq = sequence( trav_path );
          auto it = trav_seq.begin() + offset;
          end = trav_seq.begin() + label.size();
          if ( end > trav_seq.end() - seed_len + 1 ) end = trav_seq.end() - seed_len + 1;
          for ( std::size_t i = offset; it != end; ++it, ++i ) {
            if ( bv[ i ] == 1 ) continue;
            range = index.find( it, it + seed_len );
            // :TODO:Wed Aug 07 11:09:\@cartoonist:
            // Getting non-empty range doesn't necessarily mean that the locus is covered.
            // The locations corresponding to the range should also be verified.
            if ( gcsa::Range::empty( range ) ) bv[ i ] = 1;
          }
          while ( offset < label.size() && bv[ offset ] == 1 ) ++offset;
        }
        --bt_itr;
        trim_back( trav_path, *bt_itr );
      }
      for ( std::size_t i = 0; i < bv.size(); ++i ) {
        if ( bv[ i ] == 1 ) {
          ++uncovered_loci;
          bv[ i ] = 0;
        }
      }
      done_idx += label.size();
      grem::clear( trav_path );
    }
  }
  std::cout << "Done in " << Timer<>::get_lap_str( "query" ) << "." << std::endl;
  double uncovered_ratio = static_cast<double>( ::uncovered_loci ) / ::total_no;
  std::cout << "Found " << ::uncovered_loci << " uncovered loci out of " << ::total_no
            << " (" << std::setprecision( 2 ) << uncovered_ratio * 100 << "%)."
            << std::endl;
}


  inline seqan::ArgumentParser::ParseResult
parse_args( Options& options, int argc, char* argv[] )
{
  // setup ArgumentParser.
  seqan::ArgumentParser parser( release::name );
  setup_argparser( parser );
  // Embedding program's meta data and build information.
  setShortDescription( parser, release::short_desc );
  setVersion( parser, release::version );
  setDate( parser, UPDATE_DATE );
  addDescription( parser, release::desc );
  // parse command line.
  auto res = seqan::parse( parser, argc, argv );
  // only extract options if the program will continue after parse_args()
  if ( res != seqan::ArgumentParser::PARSE_OK ) {
    return res;
  }

  get_option_values( options, parser );

  return seqan::ArgumentParser::PARSE_OK;
}


  inline void
setup_argparser( seqan::ArgumentParser& parser )
{
  // positional arguments.
  std::string POSARG1 = "GRAPH_FILE";
  // add usage line.
  addUsageLine( parser, "[\\fIOPTIONS\\fP] \"\\fI" + POSARG1 + "\\fP\"" );
  // sequence file -- positional argument.
  seqan::ArgParseArgument vg_arg( seqan::ArgParseArgument::INPUT_FILE, POSARG1 );
  addArgument( parser, vg_arg );
  // GCSA2 index file -- **required** option.
  seqan::ArgParseOption gcsa_arg( "g", "gcsa", "GCSA2 index file.",
      seqan::ArgParseArgument::INPUT_FILE, "GCSA2_FILE" );
  setValidValues( gcsa_arg, gcsa::GCSA::EXTENSION );
  addOption( parser, gcsa_arg );
  setRequired( parser, "g" );
  // Seed length.
  addOption( parser, seqan::ArgParseOption( "l", "seed-len", "Seed length.",
        seqan::ArgParseArgument::INTEGER, "INT" ) );
  setRequired( parser, "l" );
  // Output file.
  seqan::ArgParseOption output_arg( "o", "output",
      "Write positions where sequences are matched.",
      seqan::ArgParseArgument::OUTPUT_FILE, "OUTPUT" );
  addOption( parser, output_arg );
  setRequired( parser, "o" );
}


  inline void
get_option_values( Options& options, seqan::ArgumentParser& parser )
{
  getArgumentValue( options.vg_filename, parser, 0 );
  getOptionValue( options.gcsa_filename, parser, "gcsa" );
  getOptionValue( options.output_filename, parser, "output" );
  getOptionValue( options.seed_len, parser, "seed-len" );
}

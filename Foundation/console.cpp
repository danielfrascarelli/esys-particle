/////////////////////////////////////////////////////////////
//                                                         //
// Copyright (c) 2003-2017 by The University of Queensland //
// Centre for Geoscience Computing                         //
// http://earth.uq.edu.au/centre-geoscience-computing      //
//                                                         //
// Primary Business: Brisbane, Queensland, Australia       //
// Licensed under the Open Software License version 3.0    //
// http://www.apache.org/licenses/LICENSE-2.0              //
//                                                         //
/////////////////////////////////////////////////////////////

#include "Foundation/console.h"

//--- IO includes ---
#include <iostream>
#include <fstream>

//--- MPI includes ---
#include <mpi.h>

using std::ofstream;
using std::ios;

Console console;

/*!
  Console default constructor. Output mode is set to unbuffered
*/ 
Console::Console() 
{
	m_initialized=false;
	m_buffered=false;
	m_bufflen=0;
	m_mute=true; 
	m_has_new_data=false;
	m_verb=-1;// prevent writing to non-initialized console
}    

/*!
  Destructor - flush any remaining buffers
*/
Console::~Console() 
{
	flush();
}

/*!
  Initialize the console. This can't be done in the constructor because
  Console is intended for use as a static object and some of the needed 
  data are (MPI-Rank, time offset) are only available after MPI_Init
  Default filename is "console.out.RANK"
*/
void Console::Initialize()
{
	int rank;
	stringstream sst;	
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get MPI rank
	m_time_offs=MPI_Wtime();
	
	sst << "console.out." << rank;
 	m_filename=sst.str();
}

/*!
  Initialize the console. This can't be done in the constructor because
  Console is intended for use as a static object and some of the needed 
  data are (MPI-Rank, time offset) are only available after MPI_Init

  \param filename the name of the output file
*/
void Console::Initialize(const string& filename)
{
	int rank;
	stringstream sst;	
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get MPI rank
	m_time_offs=MPI_Wtime();
	
	sst << filename << "." << rank;
 	m_filename=sst.str();
}


/*!
  Set base of the output file name .PID will be added where PID is
  the MPI process rank 
  
  \param filename the name of the output file
*/
void Console::SetFilename(const string& filename)
{
	int rank;
	stringstream sst;	
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get MPI rank
	
	sst << filename << "." << rank;
 	m_filename=sst.str();
}

/*!
  Set buffering mode of the console. If the buffer length parameter
  is non-zero, buffered mode is set.
  
  \param bufflen the lenght of the console buffer 
*/
void Console::SetBuffered(unsigned int bufflen)
{
	m_bufflen=bufflen;
	if(bufflen>0){
			m_buffered=true;
	} else {
			m_buffered=false;
	}
}

/*!
  set current verbosity level of the console
*/
void Console::SetVerbose(int vl)
{
	m_verb=vl;
}
	
/*!
  flush the console buffer to the output file
*/
void Console::flush() 
{
	if(m_has_new_data) {
		ofstream outfile(m_filename.c_str(),ios::app); // make sure the file is opened in append mode
		outfile << m_buffer.str();
		outfile.close();
		m_buffer.str(string(""));
		m_has_new_data=false;
	}
}

Console& Console::Critical() 
{
	if ( m_verb > 0) {
		m_mute = false ;
		double curr_time = MPI_Wtime()-m_time_offs;
		m_buffer << "[[ " << curr_time << " ]] crt - " ;
		m_has_new_data=true;
	} else {
		m_mute = true ;
	}
	return *this ;
}

Console & Console::Error() 
{
	if ( m_verb > 1) {
		m_mute = false ;
		double curr_time = MPI_Wtime()-m_time_offs;
		m_buffer << "[[ " << curr_time << " ]] err - " ;
		m_has_new_data=true;
	} else {
		m_mute = true ;
	}
	return *this ;
}

Console & Console::Warning() 
{
	if ( m_verb > 2) {
		m_mute = false ;
		double curr_time = MPI_Wtime()-m_time_offs;
		m_buffer << "[[ " << curr_time << " ]] wrn - " ;
		m_has_new_data=true;
	} else {
		m_mute = true ;
	}
	return *this ; 
}

Console & Console::Message() 
{
	if ( m_verb > 3) {
		m_mute = false ;
		double curr_time = MPI_Wtime()-m_time_offs;
		m_buffer << "[[ " << curr_time << " ]] msg - " ;
		m_has_new_data=true;
	} else {
		m_mute = true ;
	}
	return *this ; 
}

Console & Console::Info() 
{
	if ( m_verb > 4) {
		m_mute = false ;
		double curr_time = MPI_Wtime()-m_time_offs;
		m_buffer << "[[ " << curr_time << " ]] inf - " ;
		m_has_new_data=true;
	} else {
		m_mute = true ;
	}
	return *this ;   
}

Console & Console::Debug() 
{
	if ( m_verb > 5) {
		m_mute = false ;
		double curr_time = MPI_Wtime()-m_time_offs;
		m_buffer << "[[ " << curr_time << " ]] dbg - " ;
		m_has_new_data=true;
	} else {
		m_mute = true ;
	}
	return *this ;  
}

Console & Console::XDebug() 
{
	if ( m_verb > 6) {
		m_mute = false ;
		double curr_time = MPI_Wtime()-m_time_offs;
		m_buffer << "[[ " << curr_time << " ]] xbg - " ;
		m_has_new_data=true;
	} else {
		m_mute = true ;
	}
	return *this ;  
}



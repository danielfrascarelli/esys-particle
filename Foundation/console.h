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

#ifndef _CONSOLE_H_
#define _CONSOLE_H_

//--STL includes--
#include <string>
#include <sstream>

using std::string;
using std::stringstream;

/*!
  \class Console
  \brief Handle message ouput on the console 
*/
class Console 
{
protected:
	int m_verb ;           //!< current verbose level
	bool m_mute;      	//!< flag to mute console (loats output) or quiet console
	bool m_buffered;
	bool m_initialized;	//!< set if filename and so on are set
	bool m_has_new_data; 
	
	unsigned int m_bufflen;	//!< length of internal buffer

	stringstream m_buffer;	//!< buffer for output & process
	string m_filename;        //!< Output stream 

	double m_time_offs;	//!< time offset

	void flush();

public:
	Console();
	virtual ~Console() ;
    
	void Initialize();
	void Initialize(const string&);
	void SetBuffered(unsigned int); //!< set buffer length and buffered mode on/off
	void SetFilename(const string&);
	void SetVerbose(int vl=7) ; //!< set verbose level - defaults to all
	inline int GetVerbose() { return m_verb; } ;

	Console & Message() ;  //!< set verbose level of next message to "msg"
	Console & Error() ;    //!< set verbose level of next message to "err"
	Console & Warning() ;  //!< set verbose level of next message to "wrn"
	Console & Critical() ; //!< set verbose level of next message to "crt"
	Console & Info() ;     //!< set verbose level of next message to "inf"
	Console & Debug() ;    //!< set verbose level of next message to "dbg"
	Console & XDebug() ;   //!< set verbose level of next message to "xdg"
    
	template <class T> Console &  operator<<(T);

} ;

#include "console.hpp"

extern Console console;

#endif // _CONSOLE_H_

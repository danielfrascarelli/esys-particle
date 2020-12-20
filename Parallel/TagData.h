#ifndef __TAGDATA_H
#define __TAGDATA_H

//--- Project includes ---

#include "Parallel/mpivbuf.h"
#include "Foundation/console.h"

/*!
	Utility class for the packing / unpacking of interaction tag data
	to / from an mpi buffer 
*/
class CTagData 
{
	private:
		int m_tag1, m_mask1, m_tag2 , m_mask2;
	
	public:
		CTagData();
		
		void setFromMPIBuffer(CVarMPIBuffer& param_buffer);
	
		inline int Tag1() const {return m_tag1;};
		inline int Tag2() const {return m_tag2;};
		inline int Mask1() const {return m_mask1;};
		inline int Mask2() const {return m_mask2;};
		
};

#endif //__TAGDATA_H
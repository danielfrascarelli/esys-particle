//--- Project includes ---

#include "TagData.h"

/*!
	default constructor - sets tags & masks to 0
*/
CTagData::CTagData(): m_tag1(0),m_tag2(0),m_mask1(0),m_mask2(0)
{}
	
/*!
	unpack data from an MPI buffer 

	\param param_buffer the buffer
*/
void CTagData::setFromMPIBuffer(CVarMPIBuffer& param_buffer)
{
	m_tag1=param_buffer.pop_int();
	m_mask1=param_buffer.pop_int();
	m_tag2=param_buffer.pop_int();
	m_mask2=param_buffer.pop_int();
	console.XDebug() << "TagData: tag1, mask1, tag2, mask2 " 
		<< m_tag1 << " , " << m_mask1 << " , " 
		<< m_tag2 << " , " << m_mask2 << "\n";  
}
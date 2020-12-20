template <class T>
Console&  Console::operator<<(T payload) 
{
	// write payload to buffer 
	if (!m_mute){
		m_buffer << payload;
	}
    
	// if unbuffered mode or buffer length exceeded, flush buffer 
	if(!m_buffered || (m_buffer.str().length() > m_bufflen)){
		flush();
	}
	
	return *this ;
}

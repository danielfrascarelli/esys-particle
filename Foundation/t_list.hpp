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

template <class T>
void Stack<T>::Push(T* V) {
  L.InsertAtStart(V) ;
} 

template <class T>
T* Stack<T>::Pop() {
  L.First() ;
  T* t= L.Get() ;
  L.Clear() ;
  return t ;
}

template <class T>
List<T>::List()
   {
   Start = End = Current = NULL ;
   }

template <class T>
List<T>::List(const List<T> &L)
   {
   Node<T> *Lcur ;
   Start = End = Current = NULL ;
   Lcur = L.Start ;
   while (Lcur)
      {
      Append(Lcur->Val) ;
      Lcur = Lcur->Next ;
      }
   }

template <class T>
List<T>::~List()
   {
   Destroy() ;
   }

/*!
  return the number of element in the list
*/
template <class T>
int List<T>::SizeList()
   {
    int i = 0;
    First() ;
    while (!IsEnd()) {
	i++ ;
	Next() ;
    }
    return i ;
   }

/*!
  Insert an element at the begining of the list
  \param V pointer to element to be inserted
*/
template <class T>
void List<T>::InsertAtStart(T *V)
   {
   if (Start==NULL)
      {
      Start = new Node<T> ;
      Start->Next = NULL ;
      Start->Prev = NULL ;
      Current = Start ;
      End = Start ;
      Current->Val = V ;
      }
   else
      {
      Node<T> *Inter ;
      Inter = new Node<T> ;
      Inter->Prev = NULL ;
      Inter->Next = Start ;
      Inter->Val = V ;
      Start->Prev = Inter ;
      Start = Inter ;
      }
   }

template <class T>
void List<T>::Append(T *V)
   {
   if (Start==NULL)
      {
      Start = new Node<T> ;
      Start->Next = NULL ;
      Start->Prev = NULL ;
      Current = Start ;
      End = Start ;
      Current->Val = V ;
      }
   else
      {
      Node<T> *Inter ;
      Inter = new Node<T> ;
      Inter->Prev = End ;
      Inter->Next = NULL ;
      Inter->Val = V ;
      End->Next = Inter ;
      End = Inter ;
      }
   }

template <class T>
void List<T>::InsertAfter(T *V)
   {
   if ((Start==NULL)||(Current==NULL)||(Current->Next==NULL))
      {
      Append(V) ;
      }
   else
      {
      Node<T> *Inter ;
      Inter = new Node<T> ;
      Inter->Prev = Current ;
      Inter->Next = Current->Next ;
      Inter->Val = V ;
      Current->Next->Prev = Inter ;
      Current->Next = Inter ;
      }
   }

template <class T>
void List<T>::InsertBefore(T *V)
   {
   if ((Start==NULL)||(Current==NULL))
      {
      Append(V) ;
      }
   else if (Start==Current)
      {
      InsertAtStart(V) ;
      }
   else
      {
      Node<T> *Inter ;
      Inter = new Node<T> ;
      Inter->Prev = Current->Prev ;
      Inter->Next = Current ;
      Inter->Val = V ;
      Current->Prev->Next = Inter ;
      Current->Prev = Inter ;
      }
   }

template <class T>
T* List<T>::Get()
   {
#ifdef _DEBUG
     if (!Current) {
       console.Critical() << "Invalid access in List via null pointer\n" ;
     }
     if (!Current->Val) {
       console.Warning() << "List just returned a null pointer\n" ;
     }
#endif
   return Current->Val ;
   }

template <class T>
void List<T>::Put(T *V)
   {
   if (Current) Current->Val = V ;
   }

template <class T>
void List<T>::Clear()
   {
   if (Current)
      {
      if (Current==Start)
	 {
	 Start = Current->Next ;
	 if (Start) Start->Prev = NULL ;
	 else End=NULL ;
	 delete Current ;
	 Current = Start ;
	 }
      else if (Current==End)
	 {
	 End = Current->Prev ;
	 End->Next = NULL ;
	 delete Current ;
	 Current = NULL ;
	 }
      else
	 {
	 Node<T> *Inter ;

	 Inter = Current->Next ;
	 Current->Prev->Next = Current->Next ;
	 Current->Next->Prev = Current->Prev ;
	 delete Current ;
	 Current = Inter ;
	 }
      }
   }

template <class T>
void List<T>::Destroy()
   {
   Current=Start ;
   while (Current) Clear() ;
   }

template <class T>
List<T>& List<T>::operator << (T *V)
  {
  InsertAfter(V) ;
  Next() ;
  return *this ;
  }

template <class T>
List<T>& List<T>::operator >> (T *V)
  {
  if (Current) V = Current->Val ;
  Next() ;
  return *this ;
  }

template <class T>
void List<T>::Next()
  {
  if (Current) Current=Current->Next ;
  }

template <class T>
void List<T>::Prev()
  {
  if (Current && Current->Prev) Current=Current->Prev ;
  }

template <class T>
void List<T>::First()
   {
   Current = Start ;
   }

template <class T>
void List<T>::Last()
   {
   Current = End ;
   }

template <class T>
int List<T>::IsEnd()
   {
   return (Current==NULL) ;
   }

template <class T>
int List<T>::IsStart()
   {
   return (Current==Start) ;
   }

template <class T>
void List<T>::Swap() // echange current avec le suivant
   {
   Node<T> *CS, *CP, *CSS, *C ;

   C = Current ;
   CS = C->Next ;
   CP = C->Prev ;
   CSS = CS->Next ;
   if (CP) CP->Next = CS ;
   CS->Prev = CP ;
   CS->Next = C ;
   C->Prev = CS ;
   C->Next = CSS ;
   if (CSS) CSS->Prev = C ;
   if (!CP) Start = CS ;
   if (!CSS) End = C ;
   }

template <class T>
List<T> List<T>::operator + (const List<T>& L)
   {
   List<T> This=*this;
   This+=L ;
   return This ;
   }

template <class T>
List<T>& List<T>::operator += (const List<T>& L)
   {
   Node<T> *Lcur ;
   Lcur = L.Start ;
   while (Lcur)
      {
      Append(Lcur->Val) ;
      Lcur = Lcur->Next ;
      }
   return *this ;
   }

template <class T>
List<T>& List<T>::operator = (const List<T>& L)
   {
   Node<T> *Lcur ;
   Destroy() ;
   Lcur = L.Start ;
   while (Lcur)
      {
      Append(Lcur->Val) ;
      Lcur = Lcur->Next ;
      }
   return *this ;
   }

#if 0
template <class T>
void ListWS<T>::PushPos() 
{ 
  stack.Push(Current) ;
}
 
template <class T>
void ListWS<T>::PopPos() 
{ 
  Current = stack.Pop() ;
} 
#endif

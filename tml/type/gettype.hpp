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

//=== PAIR ==== 
template<typename T1,typename T2>
MPI_Datatype SGetType::operator()(const pair<T1,T2>& P)
{
  if(!tml_pair<T1,T2>::initialized){
    MPI_Aint addr_first,addr_second;
    MPI_Aint disp[2];
    MPI_Datatype type[2];
    int blocklen[2]={1,1};

    MPI_Get_address((void*)&P.first,&addr_first);
    MPI_Get_address((void*)&P.second,&addr_second);
    disp[0]=MPI_Aint(0);
    disp[1]=addr_second-addr_first;
    type[0]=GetType(P.first);
    type[1]=GetType(P.second);
    MPI_Type_create_struct(2,blocklen,disp,type,&tml_pair<T1,T2>::type);
    MPI_Type_commit(&tml_pair<T1,T2>::type);
    tml_pair<T1,T2>::initialized=true;
  }
  return tml_pair<T1,T2>::type;
}

// === TRIPLET===
template<typename T1, typename T2, typename T3>
MPI_Datatype SGetType::operator()(const esys::lsm::triplet<T1,T2,T3> &PP)
{
  if(!tml_trip<T1,T2,T3>::initialized){
    esys::lsm::triplet<T1,T2,T3> P(PP);
    MPI_Aint addr[3];
    MPI_Aint disp[3];
    MPI_Datatype type[3];
    int blocklen[3]={1,1,1};

    MPI_Get_address(static_cast<void *>(&(P.template get<0>())), &(addr[0]));
    MPI_Get_address(static_cast<void *>(&(P.template get<1>())), &(addr[1]));
    MPI_Get_address(static_cast<void *>(&(P.template get<2>())), &(addr[2]));
    disp[0]=MPI_Aint(0);
    disp[1]=addr[1]-addr[0];
    disp[2]=addr[2]-addr[0];
    type[0]=GetType(P.template get<0>());
    type[1]=GetType(P.template get<1>());
    type[2]=GetType(P.template get<2>());
    MPI_Type_create_struct(3,blocklen,disp,type,&tml_trip<T1,T2,T3>::type);
    MPI_Type_commit(&tml_trip<T1,T2,T3>::type);
    tml_trip<T1,T2,T3>::initialized = true;
  }
  return tml_trip<T1,T2,T3>::type;
}

//=== QUAD ===
template<typename T1, typename T2, typename T3, typename T4>
MPI_Datatype SGetType::operator()(const esys::lsm::quadtuple<T1,T2,T3,T4> &PP)
{
  if(!tml_quad<T1,T2,T3,T4>::initialized){
    esys::lsm::quadtuple<T1,T2,T3,T4> P(PP);
    MPI_Aint addr[4];
    MPI_Aint disp[4];
    MPI_Datatype type[4];
    int blocklen[4]={1,1,1,1};

    MPI_Get_address(static_cast<void *>(&(P.template get<0>())), &(addr[0]));
    MPI_Get_address(static_cast<void *>(&(P.template get<1>())), &(addr[1]));
    MPI_Get_address(static_cast<void *>(&(P.template get<2>())), &(addr[2]));
    MPI_Get_address(static_cast<void *>(&(P.template get<3>())), &(addr[3]));
    disp[0]=MPI_Aint(0);
    disp[1]=addr[1]-addr[0];
    disp[2]=addr[2]-addr[0];
    disp[3]=addr[3]-addr[0];
    type[0]=GetType(P.template get<0>());
    type[1]=GetType(P.template get<1>());
    type[2]=GetType(P.template get<2>());
    type[3]=GetType(P.template get<3>());
    MPI_Type_create_struct(4,blocklen,disp,type,&tml_quad<T1,T2,T3,T4>::type);
    MPI_Type_commit(&tml_quad<T1,T2,T3,T4>::type);
    tml_quad<T1,T2,T3,T4>::initialized = true;
  }
  return tml_quad<T1,T2,T3,T4>::type;
}

//=== QUINT ===
template<typename T1, typename T2, typename T3, typename T4, typename T5>
MPI_Datatype SGetType::operator()(const esys::lsm::quintuple<T1,T2,T3,T4,T5> &PP)
{
  if(!tml_quin<T1,T2,T3,T4,T5>::initialized){
    esys::lsm::quintuple<T1,T2,T3,T4,T5> P(PP);
    MPI_Aint addr[5];
    MPI_Aint disp[5];
    MPI_Datatype type[5];
    int blocklen[5]={1,1,1,1,1};

    MPI_Get_address(static_cast<void *>(&(P.template get<0>())), &(addr[0]));
    MPI_Get_address(static_cast<void *>(&(P.template get<1>())), &(addr[1]));
    MPI_Get_address(static_cast<void *>(&(P.template get<2>())), &(addr[2]));
    MPI_Get_address(static_cast<void *>(&(P.template get<3>())), &(addr[3]));
    MPI_Get_address(static_cast<void *>(&(P.template get<4>())), &(addr[4]));
    disp[0]=MPI_Aint(0);
    disp[1]=addr[1]-addr[0];
    disp[2]=addr[2]-addr[0];
    disp[3]=addr[3]-addr[0];
    disp[4]=addr[4]-addr[0];
    type[0]=GetType(P.template get<0>());
    type[1]=GetType(P.template get<1>());
    type[2]=GetType(P.template get<2>());
    type[3]=GetType(P.template get<3>());
    type[4]=GetType(P.template get<4>());
    MPI_Type_create_struct(5,blocklen,disp,type,&tml_quin<T1,T2,T3,T4,T5>::type);
    MPI_Type_commit(&tml_quin<T1,T2,T3,T4,T5>::type);
    tml_quin<T1,T2,T3,T4,T5>::initialized = true;
  }
  return tml_quin<T1,T2,T3,T4,T5>::type;
}

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

#ifndef __SUBLATTICE_CMD_H
#define __SUBLATTICE_CMD_H

const int CMD_FINISH=1024;
const int CMD_PRINT=1025;
const int CMD_CALC=1026;
const int CMD_XCHANGE=1027;
const int CMD_NSEARCH=1028;
const int CMD_UPDATE=1029;
const int CMD_CHECKNEIGHBORS=1030;
const int CMD_PMOVE=1031;
const int CMD_WALL=1033;
const int CMD_WMOVE=1034;
const int CMD_COUNT=1035;
const int CMD_INITLATTICE=1036;
const int CMD_INIT2DTRILOCAL=1037;
const int CMD_ADDEWALLIG=1038;
const int CMD_ADDBONDEDIG=1039;
const int CMD_MAKELATTICE=1041;
const int CMD_ADDPIG=1042;
const int CMD_EXIG=1043;
const int CMD_INIT2DRANDOM=1044;
const int CMD_ADDDAMP=1045;
const int CMD_PVEL=1046;
const int CMD_INIT3DTRILOCAL=1047;
const int CMD_ADDSIG=1048;
const int CMD_INITCOMPLEX=1049;
const int CMD_ADDBWALLIG=1050;
const int CMD_RECEIVEPARTICLES=1051;
const int CMD_INITLATTICECIRC=1052;
const int CMD_PSETND=1053;
const int CMD_SAVECHECKPOINT=1054;
const int CMD_GETSOFTBONDEDWALL=1055;
const int CMD_WFORCE=1056;
const int CMD_PTAG=1057;
const int CMD_GETVISCWALL=1058;
const int CMD_WVEL=1059;
const int CMD_PERFORMTIMING = 1060;
const int CMD_SAVETIMINGDATA = 1061;
const int CMD_ADDROTBONDEDIG=1062;
const int CMD_GETNUMPARTICLES=1063;
const int CMD_ADDROTDAMP=1064;
const int CMD_ADDTRIMESHIG = 1065;
const int CMD_MOVENODE = 1066;
const int CMD_MOVETAGGEDNODES = 1067;
const int CMD_DO2DCALCULATIONS = 1068;
const int CMD_ADDTRIMESH = 1069;
const int CMD_ADDBONDEDTRIMESHIG = 1071;
const int CMD_PSETNR=1072;
const int CMD_ADDMESH2D = 1073;
const int CMD_ADDBONDEDMESH2DIG = 1074;
const int CMD_FINDNEARESTPARTICLE = 1075;
const int CMD_GETPARTICLEPOSN = 1076;
const int CMD_GETIDPARTICLEDATA = 1077;
const int CMD_IDPARTICLEMOVE = 1078;
const int CMD_PANGVEL = 1079;
const int CMD_ADDWALL = 1080;
const int CMD_ADDVWALLIG = 1081;
const int CMD_ADDBBWALLIG =1082;
const int CMD_SETTIMESTEPSIZE =1083;
const int CMD_ADDSHORTBONDEDIG=1084;
const int CMD_ADDCAPPEDBONDEDIG=1085;
const int CMD_ADDTAGPIG=1086;
const int CMD_ADDMESH2DIG=1087;
const int CMD_PDENS=1088;
const int CMD_REMOVEIG=1089;
const int CMD_SAVESNAPSHOT=1090;
const int CMD_LOADCHECKPOINT=1091;
const int CMD_SAVECHECKPOINTWTM=1092;
const int CMD_TRANSLATEMESHBY=1093;
const int CMD_PMOVETAGGEDBY=1094;
const int CMD_GETWALLPOS=1095;
const int CMD_RECEIVECONNECTIONS=1096;
const int CMD_PTVEL = 1097;
const int CMD_SETVERBOSE = 1099;
const int CMD_PSETNT = 1100;
const int CMD_WNORM = 1101;
const int CMD_GETWALLFORCE=1102;
const int CMD_ADDROTTHERMBONDEDIG=1103;
const int CMD_GETBONDGROUPIDPAIRS=1104;
const int CMD_SETVERBOSITY=1105;
const int CMD_ADDTAGGEDEWALLIG=1106;
const int CMD_ADDSPHEREBODY=1107;
const int CMD_ADDESPHEREBODYIG=1108;
const int CMD_GETSPHEREBODYPOS=1109;
const int CMD_GETSPHEREBODYFORCE=1110;
const int CMD_SPHEREBODYMOVE=1111;
const int CMD_INITCONSOLE=1112;
const int CMD_SETCONSOLEFNAME=1113;
const int CMD_SETCONSOLEBUFF=1114;
const int CMD_SETINTERACTIONPARAMS=1115;

/****fluid contents: begin****/
const int CMD_ADDFLUID_VEC3=1116;
const int CMD_RECV_PRESSURE=1117;
const int CMD_SOLVE=1118;
const int CMD_UPDATE_FLUID=1119;
const int CMD_EXCHANGE_CELLS=1120;
const int CMD_ADDFLUID=1121;
/****fluid contents: end****/

const int CMD_ADDBRITTLEBEAMSCIG=1122;
const int CMD_ADDBRITTLEBEAMDZCIG=1123;
const int CMD_APPLYPRESSURRETOMESH=1124;
const int CMD_ROTATEMESHBY=1125;
const int CMD_RPROT=1126;
const int CMD_CHGRADIUS=1127;

// ---------------------------------

const int CMD_ADD_SPF=2051;
const int CMD_ADD_VPF=2052;
const int CMD_SEND_FIELDS=2053;
const int CMD_ADD_SIF=2054;
const int CMD_ADD_VIF=2055;
const int CMD_ADD_VTF=2056;
const int CMD_ADD_STF=2057;
const int CMD_ADD_VWF=2058;
const int CMD_ADD_HIF=2059;
/****fluid contents: begin****/
const int CMD_ADD_SFF=2060;
const int CMD_ADD_VFF=2061;
const int CMD_ADD_SFIF=2062;
const int CMD_ADD_VFIF=2063;
const int CMD_SEND_COEFFI=2064;
/****fluid contents: end****/


// -------------------------------

const int CMD_GETMESHNODEREF=4096;
const int CMD_GETMESHFACEREF=4097;
const int CMD_GETMESH2DSTRESS=4098;
const int CMD_GETTRIMESHFORCE=4099;



#endif //__SUBLATTICE_CMD_H

/*=====================================================================================*/
/*! \file		FEMElementTri.cpp
	\author		peterkau
	\brief		Implementation of class FEMElementTri
 */
/*=====================================================================================*/

#include "SimpleFEMDefs.h"
#include "FEMElementTri.h"
#include "FEMMesh.h"
#include <cmath>

static const bool useGeom = true;

// TASK 3
void FEMElementTri::Assemble(FEMMesh *pMesh) const
{
  Vector2 v0 = pMesh->GetNodePosition(m_nodes[0]); 
  Vector2 v1 = pMesh->GetNodePosition(m_nodes[1]); 
  Vector2 v2 = pMesh->GetNodePosition(m_nodes[2]); 
  Vector2 e01,e12;
  e01 = v1-v0;
  e12 = v2-v1;
  double area = abs(e01.x()*(-e12).y() - e01.y()*(-e12).x()) / 2;

  Vector2 d_N_d_x[3];
  for(int i=0; i<3; i++){
    if(useGeom)
      computeSingleBasisDerivGlobalGeom(i, d_N_d_x[i], pMesh);
    else
      computeSingleBasisDerivGlobalLES(i, d_N_d_x[i], pMesh);
  }
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      int iGlobal = GetGlobalNodeForElementNode(i);
      int jGlobal = GetGlobalNodeForElementNode(j);
      if(iGlobal>=jGlobal){
        pMesh->AddToStiffnessMatrix(iGlobal, jGlobal, area*(d_N_d_x[i]|d_N_d_x[j]));
      }
    }
  }
}

// TASK 2
void FEMElementTri::computeSingleBasisDerivGlobalGeom(size_t nodeId, Vector2 &basisDerivGlobal, const FEMMesh *pMesh) const
{
  Vector2 v0 = pMesh->GetNodePosition(m_nodes[0]); 
  Vector2 v1 = pMesh->GetNodePosition(m_nodes[1]); 
  Vector2 v2 = pMesh->GetNodePosition(m_nodes[2]); 

  Vector2 e,e01,e12,e20;
  e01 = v1-v0;
  e12 = v2-v1;
  e20 = v0-v2;

  double area = abs(e01.x()*(-e12).y() - e01.y()*(-e12).x()) / 2;

  double a, b, height;
  switch(nodeId){
    case 0:
      e = e12;
      break;
    case 1:
      e = e20;
      break;
    case 2:
      e = e01;
      break;
  }
  a = -e.y();
  b = e.x();
  height = 2 * area / e.length();

  basisDerivGlobal = Vector2(a,b).normalized() / height;
}

// TASK 1
// nodeId is local(0~2) 
void FEMElementTri::computeSingleBasisDerivGlobalLES(size_t nodeId, Vector2 &basisDerivGlobal, const FEMMesh *pMesh) const
{
  Matrix3x3 basisMat;
  for(int i=0; i<3; i++){
    Vector2 nodePos = pMesh->GetNodePosition(m_nodes[i]);
    basisMat(i,0) = nodePos[0];
    basisMat(i,1) = nodePos[1];
    basisMat(i,2) = 1;
  }
  Vector3 delta(0.0, 0.0, 0.0);
  delta[nodeId] = 1.0;
  Vector3 abc =  basisMat.inverse() * delta;
  basisDerivGlobal = Vector2(abc[0], abc[1]);
	
}

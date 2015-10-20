/*=====================================================================================*/
/*! \file		FEMElementTri.cpp
	\author		peterkau
	\brief		Implementation of class FEMElementTri
 */
/*=====================================================================================*/

#include "SimpleFEMDefs.h"
#include "FEMElementTri.h"
#include "FEMMesh.h"
#include "cmath"
// TASK 3
void FEMElementTri::Assemble(FEMMesh *pMesh) const
{
    double value =0.,iGlobal=0.,jGlobal=0.;
    Vector2 dN_i,dN_j;
    double A=0.;
    Vector2 l1=pMesh->GetNodePosition(GetGlobalNodeForElementNode(0))-pMesh->GetNodePosition(GetGlobalNodeForElementNode(1));
    Vector2 l2=pMesh->GetNodePosition(GetGlobalNodeForElementNode(0))-pMesh->GetNodePosition(GetGlobalNodeForElementNode(2));
    
    if (l1.length()>l2.length())
    {
        A=l2.length()*l2.length()/2;
    }else
    {
        A=l1.length()*l1.length()/2;
    }
    
    for(int i=0; i<=2;++i)
    {
        computeSingleBasisDerivGlobalLES(i,dN_i,pMesh);
        iGlobal = GetGlobalNodeForElementNode(i);
        
        for (int j=0;j<=2; ++j)
        {
            jGlobal = GetGlobalNodeForElementNode(j);
            computeSingleBasisDerivGlobalLES(j,dN_j,pMesh);
            
            if(iGlobal>=jGlobal)
            {
                value = A*(dN_i.x()*dN_j.x() + dN_i.y()*dN_j.y());
                pMesh->AddToStiffnessMatrix(iGlobal, jGlobal, value);
            }
        }
    }
}

// TASK 2
void FEMElementTri::computeSingleBasisDerivGlobalGeom(size_t nodeId, Vector2 &basisDerivGlobal, const FEMMesh *pMesh) const
{
    double x=0.,y=0.,h=0.;
    Vector2 N0 = pMesh->GetNodePosition(GetGlobalNodeForElementNode(0));
    Vector2 N1 = pMesh->GetNodePosition(GetGlobalNodeForElementNode(1));
    Vector2 N2 = pMesh->GetNodePosition(GetGlobalNodeForElementNode(2));
    
    Vector2 E01= N1-N0;
    Vector2 E12= N2-N1;
    Vector2 E20= N0-N2;
    
    if(E01.length()<=E12.length())
    {
        h=sqrt(2*E01.length()*E01.length());
    }else {
        h=sqrt(2*E12.length()*E12.length());
    }
    
    switch (nodeId){
            case 0:
            x=E12.y();
            y=-E12.x();
            break;
            
            case 1:
            x=E20.y();
            y=-E20.x();
            break;
            
            case 2:
            x=E01.y();
            y=-E01.x();
            break;
    }
    
    
    basisDerivGlobal = Vector2(x,y).normalized()*h;
    
}

// TASK 1
void FEMElementTri::computeSingleBasisDerivGlobalLES(size_t nodeId, Vector2 &basisDerivGlobal, const FEMMesh *pMesh) const
{
    
    Vector3 x(pMesh->GetNodePosition(GetGlobalNodeForElementNode(0)).x(),pMesh->GetNodePosition(GetGlobalNodeForElementNode(1)).x(),pMesh->GetNodePosition(GetGlobalNodeForElementNode(2)).x());
    Vector3 y(pMesh->GetNodePosition(GetGlobalNodeForElementNode(0)).y(),pMesh->GetNodePosition(GetGlobalNodeForElementNode(1)).y(),pMesh->GetNodePosition(GetGlobalNodeForElementNode(2)).y());
    Matrix3x3 M;
    
    M(0,2)=M(1,2)=M(2,2)=1.;
    for (int i=0; i<=2; ++i) {
        M(i,0)=x[i];
        M(i,1)=y[i];
    }
    
    Vector3 n_x(0.,0.,0.);
    n_x[nodeId]=1.;
    
    Vector3 abc = M.inverse()*n_x;
    
    basisDerivGlobal = Vector2(abc.x(),abc.y());
}

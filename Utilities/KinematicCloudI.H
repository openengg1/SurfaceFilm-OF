/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fvmSup.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::KinematicCloud<CloudType>&
Foam::KinematicCloud<CloudType>::cloudCopy() const
{
    return cloudCopyPtr_();
}


template<class CloudType>
inline const Foam::fvMesh& Foam::KinematicCloud<CloudType>::mesh() const
{
    return mesh_;
}


template<class CloudType>
inline const Foam::IOdictionary&
Foam::KinematicCloud<CloudType>::particleProperties() const
{
    return particleProperties_;
}


template<class CloudType>
inline const Foam::IOdictionary&
Foam::KinematicCloud<CloudType>::outputProperties() const
{
    return outputProperties_;
}


template<class CloudType>
inline Foam::IOdictionary& Foam::KinematicCloud<CloudType>::outputProperties()
{
    return outputProperties_;
}


template<class CloudType>
inline const Foam::cloudSolution&
Foam::KinematicCloud<CloudType>::solution() const
{
    return solution_;
}


template<class CloudType>
inline Foam::cloudSolution& Foam::KinematicCloud<CloudType>::solution()
{
    return solution_;
}


template<class CloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::KinematicCloud<CloudType>::constProps() const
{
    return constProps_;
}


template<class CloudType>
inline const Foam::dictionary&
Foam::KinematicCloud<CloudType>::subModelProperties() const
{
    return subModelProperties_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::KinematicCloud<CloudType>::rho() const
{
    return rho_;
}



template<class CloudType>
inline const Foam::volVectorField& Foam::KinematicCloud<CloudType>::U() const
{
    return U_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::KinematicCloud<CloudType>::mu() const
{
    return mu_;
}

template<class CloudType>
inline const Foam::volVectorField& Foam::KinematicCloud<CloudType>::vor() const
{
    return vor_;
}


template<class CloudType>
inline const Foam::dimensionedVector& Foam::KinematicCloud<CloudType>::g() const
{
    return g_;
}


template<class CloudType>
inline Foam::scalar Foam::KinematicCloud<CloudType>::pAmbient() const
{
    return pAmbient_;
}


template<class CloudType>
inline Foam::scalar& Foam::KinematicCloud<CloudType>::pAmbient()
{
    return pAmbient_;
}


template<class CloudType>
//inline const typename CloudType::parcelType::forceType&
inline const typename Foam::KinematicCloud<CloudType>::forceType&
Foam::KinematicCloud<CloudType>::forces() const
{
    return forces_;
}


template<class CloudType>
inline typename Foam::KinematicCloud<CloudType>::functionType&
Foam::KinematicCloud<CloudType>::functions()
{
    return functions_;
}


template<class CloudType>
inline const Foam::InjectionModelList<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::injectors() const
{
    return injectors_;
}


template<class CloudType>
inline Foam::InjectionModelList<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::injectors()
{
    return injectors_;
}


template<class CloudType>
inline const Foam::DispersionModel<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::dispersion() const
{
    return dispersionModel_;
}


template<class CloudType>
inline Foam::DispersionModel<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::dispersion()
{
    return dispersionModel_();
}


template<class CloudType>
inline const Foam::PatchInteractionModel<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::patchInteraction() const
{
    return patchInteractionModel_;
}


template<class CloudType>
inline Foam::PatchInteractionModel<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::patchInteraction()
{
    return patchInteractionModel_();
}


template<class CloudType>
inline const Foam::SurfaceFilmModel<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::surfaceFilm() const
{
    return surfaceFilmModel_();
}


template<class CloudType>
inline Foam::SurfaceFilmModel<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::surfaceFilm()
{
    return surfaceFilmModel_();
}


template<class CloudType>
inline const Foam::vectorIntegrationScheme&
Foam::KinematicCloud<CloudType>::UIntegrator() const
{
    return UIntegrator_;
}

template<class CloudType>
inline const Foam::vectorIntegrationScheme&
Foam::KinematicCloud<CloudType>::OmegaIntegrator() const
{
    return OmegaIntegrator_;
}

template<class CloudType>
inline Foam::label Foam::KinematicCloud<CloudType>::nParcels() const
{
    return this->size();
}


template<class CloudType>
inline Foam::scalar Foam::KinematicCloud<CloudType>::massInSystem() const
{
    scalar sysMass = 0.0;
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
         const parcelType& p = iter();
         sysMass += p.nParticle()*p.mass();
    }

    return sysMass;
}


template<class CloudType>
inline Foam::vector
Foam::KinematicCloud<CloudType>::linearMomentumOfSystem() const
{
    vector linearMomentum(vector::zero);

    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();

        linearMomentum += p.nParticle()*p.mass()*p.U();
    }

    return linearMomentum;
}


template<class CloudType>
inline Foam::scalar
Foam::KinematicCloud<CloudType>::linearKineticEnergyOfSystem() const
{
    scalar linearKineticEnergy = 0.0;

    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();

        linearKineticEnergy += p.nParticle()*0.5*p.mass()*(p.U() & p.U());
    }

    return linearKineticEnergy;
}


template<class CloudType>
inline Foam::scalar
Foam::KinematicCloud<CloudType>::rotationalKineticEnergyOfSystem() const
{
    scalar rotationalKineticEnergy = 0.0;

    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();

        rotationalKineticEnergy +=
            p.nParticle()*0.5*p.momentOfInertia()*(p.omega() & p.omega());
    }

    return rotationalKineticEnergy;
}


template<class CloudType>
inline Foam::scalar Foam::KinematicCloud<CloudType>::Dij
(
    const label i,
    const label j
) const
{
    scalar si = 0.0;
    scalar sj = 0.0;
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        si += p.nParticle()*pow(p.d(), i);
        sj += p.nParticle()*pow(p.d(), j);
    }

    reduce(si, sumOp<scalar>());
    reduce(sj, sumOp<scalar>());
    sj = max(sj, VSMALL);

    return si/sj;
}


template<class CloudType>
inline Foam::scalar Foam::KinematicCloud<CloudType>::Dmax() const
{
    scalar d = -GREAT;
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        d = max(d, p.d());
    }

    reduce(d, maxOp<scalar>());

    return max(0.0, d);
}


// Rohit Mishra (10/01/2018)
// Code to find the impigement diameter
template<class CloudType>
inline Foam::scalar Foam::KinematicCloud<CloudType>::diameterImpig
(
const scalar radX
) const
{
vector zAxis(0,0,1);
 vector yAxis(0,1,0);
    vector xAxis(1,0,0);

scalar count = 0;
 scalar si = 0.0;
    scalar sj = 0.0;
forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        if ((p.typeId() < 1) && ( mag(p.position() & xAxis) > radX) && (mag(p.position() & xAxis) < (radX+1e-3)) && (mag(p.position() & yAxis) < 10e-3) && (mag(p.position() & zAxis) < (1e-3)))
        {    
        si += p.nParticle()*pow(p.d(), 3);
        sj += p.nParticle()*pow(p.d(), 2);
    
             count = count + 1;
          }
else
{}

    }


if(count==0)
{return 0;}

else
{


    reduce(si, sumOp<scalar>());
    reduce(sj, sumOp<scalar>());
    sj = max(sj, VSMALL);

    return si/sj;
}

}


template<class CloudType>
inline Foam::scalar Foam::KinematicCloud<CloudType>::diameterSecondImpig
(
const scalar radX
) const
{
vector zAxis(0,0,1);
 vector yAxis(0,1,0);
    vector xAxis(1,0,0);
vector pImp(10e-3,0,0);
scalar impX = mag(pImp & xAxis);
scalar count = 0;
 scalar si = 0.0;
    scalar sj = 0.0;
forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        if ((p.typeId() > 0) && ( mag(p.position() & zAxis) <(5e-3))  && (mag(p.position() & yAxis) > 5e-3) && (mag(p.position() & yAxis) < 10e-3) && (mag(p.position() & xAxis) > (impX+radX))  && (mag(p.position() & xAxis) < (impX + radX + 1e-2)))
       {
        si += p.nParticle()*pow(p.d(), 3);
        sj += p.nParticle()*pow(p.d(), 2);

             count = count + 1;
          }
else
{}

    }


if(count==0)
{return 0;}

else
{


    reduce(si, sumOp<scalar>());
    reduce(sj, sumOp<scalar>());
    sj = max(sj, VSMALL);

    return si/sj;
}

}

template<class CloudType>
inline Foam::scalar Foam::KinematicCloud<CloudType>::diameterSecondImpig_center
(
const scalar radX
) const
{
vector zAxis(0,0,1);
 vector yAxis(0,1,0);
    vector xAxis(1,0,0);
vector pImp(10e-3,0,0);
scalar impX = mag(pImp & xAxis);
scalar count = 0;
 scalar si = 0.0;
    scalar sj = 0.0;
forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        if ((p.typeId() > 0) && (mag(p.position() & zAxis) < (2.5e-3)) && (mag(p.position() & yAxis) < 5e-3) && (mag(p.position() & xAxis) > (impX+radX))  && (mag(p.position() & xAxis) < (impX + radX + 5e-3)))
       {
        si += p.nParticle()*pow(p.d(), 3);
        sj += p.nParticle()*pow(p.d(), 2);

             count = count + 1;
          }
else
{}

    }


if(count==0)
{return 0;}

else
{


    reduce(si, sumOp<scalar>());
    reduce(sj, sumOp<scalar>());
    sj = max(sj, VSMALL);

    return si/sj;
}

}






// Rohit Mishra (10/01/2018)
// Code to find the impigement velocity
template<class CloudType>
inline Foam::scalar Foam::KinematicCloud<CloudType>::velocityImpig
(
const scalar radX
) const
{

 vector yAxis(0,1,0);
    vector xAxis(1,0,0);
    vector zAxis(0,0,1);
scalar count = 0;
scalar vTotal = 0;
forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
       if ((p.typeId() < 1) && ( mag(p.position() & xAxis) > radX) && (mag(p.position() & xAxis) < (radX+1e-3)) && (mag(p.position() & yAxis) < 10e-3) && (mag(p.position() & zAxis) < (1e-3)))
         {
    Info <<endl<< "magpU"<<mag(p.U())<<endl<<"Position:"<<p.position()<<endl;        
   vTotal = vTotal + mag(p.U());
             count = count + 1;
          }

else
{}



    }


if(count==0)
{return 0;}

else
{return vTotal/count;}

}


// Code to find the impingement mass
template<class CloudType>
inline Foam::scalar Foam::KinematicCloud<CloudType>::massImpig() const
{

 vector yAxis(0,1,0);
    vector xAxis(1,0,0);
    vector zAxis(0,0,1);
scalar count = 0;
scalar vTotal = 0;
forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
       if ((p.typeId() < 1) && ( (sqr(mag(p.position() & xAxis)) + sqr(mag(p.position() & zAxis))) < 0.1e-3) && ( mag(p.position() & yAxis) < 5e-3) )
         {
vTotal = vTotal + p.mass()*p.nParticle();

    Info <<endl<< "Mass impinged:"<<vTotal<<endl;
  
             count = count + 1;
          }

else
{}



    }


if(count==0)
{return 0;}

else
{return vTotal;}

}




// Rohit Mishra (07/28/2018)
//Finding different penetration lengths of splashed parcels
//Axial
template<class CloudType>
inline Foam::scalar Foam::KinematicCloud<CloudType>::zpenetrationSplashed
(
    const scalar fraction
) const
{

if ((fraction < 0) || (fraction > 1))
   {
        FatalErrorIn
        (
            "inline Foam::scalar Foam::KinematicCloud<CloudType>::penetrationSplashed"
            "("
                "const scalar"
            ") const"
        )
            << "fraction should be in the range 0 < fraction < 1"
            << exit(FatalError);
    }

    scalar distance = 0.0;

    const label nParcel = this->size();
    globalIndex globalParcels(nParcel);
    const label nParcelSum = globalParcels.size();
    vector yAxis(0,1,0);
    vector xAxis(1,0,0);
    vector zAxis(0,0,1);
 
    if (nParcelSum == 0)
    {
        return distance;
    }

    scalar mSum = 0.0;
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {

        const parcelType& p = iter();
       if (p.typeId() > 0)
         {
              scalar m = p.nParticle()*p.mass();
              mSum +=m;
          }
    }
Info <<endl<< "Total splashed mass=" << mSum<<endl;

scalar targetMass = fraction*mSum;
scalar startPos = 1e-8;
scalar massTotal = 0;
vector pImp(0,35e-3,-6.9589e-3);

while (massTotal < targetMass)
{
    massTotal = 0.0;
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {

        const parcelType& p = iter();   //(p.U() & yAxis) > 0)
         if (  (p.typeId()>0)    &&   (mag((p.position() & zAxis) - (pImp & zAxis)) < startPos)  )
          {

             scalar m = p.nParticle()*p.mass();
             massTotal += m;
           // Info<<endl<<"p.U & yAxis="<<(p.U() & yAxis)<<endl<<"p.position() & zAxis=" <<(p.position() & zAxis)<<endl;
          }
    }

startPos+=1e-5;
}
distance = startPos;
return distance;
}

//Height
template<class CloudType>
inline Foam::scalar Foam::KinematicCloud<CloudType>::ypenetrationSplashed
(
    const scalar fraction
) const
{
if ((fraction < 0) || (fraction > 1))
    {
        FatalErrorIn
        (
            "inline Foam::scalar Foam::KinematicCloud<CloudType>::penetrationSplashed"
            "("
                "const scalar"
            ") const"
        )
            << "fraction should be in the range 0 < fraction < 1"
            << exit(FatalError);
    }

    scalar distance = 0.0;

    const label nParcel = this->size();
    globalIndex globalParcels(nParcel);
    const label nParcelSum = globalParcels.size();
    vector yAxis(0,1,0);
    vector xAxis(1,0,0);
    vector zAxis(0,0,1);

    if (nParcelSum == 0)
    {
        return distance;
    }

    scalar mSum = 0.0;
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {

        const parcelType& p = iter();
       if (p.typeId() > 0)
          {
              scalar m = p.nParticle()*p.mass();
              mSum +=m;
          }
    }
//Info <<endl<< "Total splashed mass=" << mSum<<endl;

scalar targetMass = fraction*mSum;
scalar startPos = 1e-8;
scalar massTotal = 0;
vector pImp(0,35e-3,-6.9589e-3);

while (massTotal < targetMass)
{


    massTotal = 0.0;
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {

        const parcelType& p = iter();
         if (  (p.typeId()>0)    &&   (mag((p.position() & yAxis) - (pImp & yAxis)) < startPos)  )
          {

             scalar m = p.nParticle()*p.mass();
             massTotal += m;
           //Info<<endl<<"p.U & yAxis="<<(p.U() & yAxis)<<endl<<"p.position() & yAxis=" <<(p.position() & yAxis)<<endl;
          }
    }

startPos+=1e-5;

}
distance = startPos;
return distance;
}



//Radial
template<class CloudType>
inline Foam::scalar Foam::KinematicCloud<CloudType>::xpenetrationSplashed
(
    const scalar fraction
) const
{
if ((fraction < 0) || (fraction > 1))
    {
        FatalErrorIn
        (   
            "inline Foam::scalar Foam::KinematicCloud<CloudType>::penetrationSplashed"
            "(" 
                "const scalar"
            ") const"
        )   
            << "fraction should be in the range 0 < fraction < 1"
            << exit(FatalError);
    }
    
    scalar distance = 0.0;
    
    const label nParcel = this->size();
    globalIndex globalParcels(nParcel);
    const label nParcelSum = globalParcels.size();
    vector yAxis(0,1,0);
    vector xAxis(1,0,0);

    if (nParcelSum == 0)
    {   
        return distance;
    }

    scalar mSum = 0.0;
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
      
        const parcelType& p = iter();
       if (p.typeId()> 0 )
          {      
              scalar m = p.nParticle()*p.mass();
              mSum +=m;
          }
    }  
//Info <<endl<< "Total splashed mass=" << mSum<<endl;

scalar targetMass = fraction*mSum;
scalar startPos = 1e-8;
scalar massTotal = 0;
vector pImp(0,35e-3,-6.9589e-3);

while (massTotal < targetMass)
{


    massTotal = 0.0;
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {

       const parcelType& p = iter();
         if (  (p.typeId()>0)    &&   (mag((p.position() & xAxis) - (pImp & xAxis)) < startPos) )
          {
             
             scalar m = p.nParticle()*p.mass();
             massTotal += m;
         //   Info<<endl<<"p.U & yAxis="<<(p.U() & yAxis)<<endl<<"p.position() & xAxis=" <<(p.position() & xAxis)<<endl;
          }
    }

startPos+=1e-5;

}
distance = startPos;
return distance;
}



template<class CloudType>
inline Foam::scalar Foam::KinematicCloud<CloudType>::penetration
(
    const scalar fraction
) const
{
    if ((fraction < 0) || (fraction > 1))
    {
        FatalErrorIn
        (
            "inline Foam::scalar Foam::KinematicCloud<CloudType>::penetration"
            "("
                "const scalar"
            ") const"
        )
            << "fraction should be in the range 0 < fraction < 1"
            << exit(FatalError);
    }

    scalar distance = 0.0;

    const label nParcel = this->size();
    globalIndex globalParcels(nParcel);
    const label nParcelSum = globalParcels.size();

    if (nParcelSum == 0)
    {
        return distance;
    }

    // lists of parcels mass and distance from initial injection point
    List<List<scalar> > procMass(Pstream::nProcs());
    List<List<scalar> > procDist(Pstream::nProcs());

    List<scalar>& mass = procMass[Pstream::myProcNo()];
    List<scalar>& dist = procDist[Pstream::myProcNo()];

    mass.setSize(nParcel);
    dist.setSize(nParcel);

    label i = 0;
    scalar mSum = 0.0;
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        scalar m = p.nParticle()*p.mass();
        scalar d = mag(p.position() - p.position0());
        mSum += m;

        mass[i] = m;
        dist[i] = d;

        i++;
    }

    // calculate total mass across all processors
    reduce(mSum, sumOp<scalar>());
    Pstream::gatherList(procMass);
    Pstream::gatherList(procDist);

    if (Pstream::master())
    {
        // flatten the mass lists
        List<scalar> allMass(nParcelSum, 0.0);
        SortableList<scalar> allDist(nParcelSum, 0.0);
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            SubList<scalar>
            (
                allMass,
                globalParcels.localSize(procI),
                globalParcels.offset(procI)
            ).assign(procMass[procI]);

            // flatten the distance list
            SubList<scalar>
            (
                allDist,
                globalParcels.localSize(procI),
                globalParcels.offset(procI)
            ).assign(procDist[procI]);
        }

        // sort allDist distances into ascending order
        // note: allMass masses are left unsorted
        allDist.sort();

        if (nParcelSum > 1)
        {
            const scalar mLimit = fraction*mSum;
            const labelList& indices = allDist.indices();

            if (mLimit > (mSum - allMass[indices.last()]))
            {
                distance = allDist.last();
            }
            else
            {
                // assuming that 'fraction' is generally closer to 1 than 0,
                // loop through in reverse distance order
                const scalar mThreshold = (1.0 - fraction)*mSum;
                scalar mCurrent = 0.0;
                label i0 = 0;

                forAllReverse(indices, i)
                {
                    label indI = indices[i];

                    mCurrent += allMass[indI];

                    if (mCurrent > mThreshold)
                    {
                        i0 = i;
                        break;
                    }
                }

                if (i0 == indices.size() - 1)
                {
                    distance = allDist.last();
                }
                else
                {
                    // linearly interpolate to determine distance
                    scalar alpha = (mCurrent - mThreshold)/allMass[indices[i0]];
                    distance =
                        allDist[i0] + alpha*(allDist[i0+1] - allDist[i0]);
                }
            }
        }
        else
        {
            distance = allDist.first();
        }
    }

    Pstream::scatter(distance);

    return distance;
}


template<class CloudType>
inline Foam::cachedRandom& Foam::KinematicCloud<CloudType>::rndGen()
{
    return rndGen_;
}


template<class CloudType>
inline Foam::List<Foam::DynamicList<typename CloudType::particleType*> >&
Foam::KinematicCloud<CloudType>::cellOccupancy()
{
    if (cellOccupancyPtr_.empty())
    {
        buildCellOccupancy();
    }

    return cellOccupancyPtr_();
}


template<class CloudType>
inline Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::KinematicCloud<CloudType>::UTrans()
{
    return UTrans_();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::KinematicCloud<CloudType>::UTrans() const
{
    return UTrans_();
}

template<class CloudType>
inline Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::KinematicCloud<CloudType>::UTransDrag()
{
    return UTransDrag_();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::KinematicCloud<CloudType>::UTransDrag() const
{
    return UTransDrag_();
}

template<class CloudType>
inline Foam::DimensionedField<Foam::tensor, Foam::volMesh>&
Foam::KinematicCloud<CloudType>::DiTrans()
{
    return DiTrans_();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::tensor, Foam::volMesh>&
Foam::KinematicCloud<CloudType>::DiTrans() const
{
    return DiTrans_();
}

template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::KinematicCloud<CloudType>::ksgsTrans()
{
    return ksgsTrans_();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::KinematicCloud<CloudType>::ksgsTrans() const
{
    return ksgsTrans_();
}




template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::KinematicCloud<CloudType>::UCoeff()
{
    return UCoeff_();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::KinematicCloud<CloudType>::UCoeff() const
{
    return UCoeff_();
}


template<class CloudType>
inline Foam::tmp<Foam::fvVectorMatrix>
Foam::KinematicCloud<CloudType>::SU(volVectorField& U) const
{
    if (debug)
    {
        Info<< "UTrans min/max = " << min(UTrans()).value() << ", "
            << max(UTrans()).value() << nl
            << "UCoeff min/max = " << min(UCoeff()).value() << ", "
            << max(UCoeff()).value() << endl;
    }

    if (solution_.coupled())
    {
        if (solution_.semiImplicit("U"))
        {
            const DimensionedField<scalar, volMesh>
                Vdt(mesh_.V()*this->db().time().deltaT());

            return UTrans()/Vdt - fvm::Sp(UCoeff()/Vdt, U) + UCoeff()/Vdt*U;
        }
        else
        {
            tmp<fvVectorMatrix> tfvm(new fvVectorMatrix(U, dimForce));
            fvVectorMatrix& fvm = tfvm();

            fvm.source() = -UTrans()/(this->db().time().deltaT());

            return tfvm;
        }
    }

    return tmp<fvVectorMatrix>(new fvVectorMatrix(U, dimForce));
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicCloud<CloudType>::theta() const
{
    tmp<volScalarField> ttheta
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":theta",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    volScalarField& theta = ttheta();
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();

        theta[cellI] += p.nParticle()*p.volume();
    }

    theta.internalField() /= mesh_.V();
    theta.correctBoundaryConditions();

    return ttheta;
}

template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicCloud<CloudType>::Vl() const
{
    tmp<volScalarField> tVl
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "Vl",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(0,3,0,0,0,0,0), scalar(0.0)),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    volScalarField& Vl = tVl();
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();

        Vl[cellI] += p.nParticle()*p.volume();
    }

    Vl.correctBoundaryConditions();

    return tVl;
}



template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicCloud<CloudType>::alpha() const
{
    tmp<volScalarField> talpha
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":alpha",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    scalarField& alpha = talpha().internalField();
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();

        alpha[cellI] += p.nParticle()*p.mass();
    }

    alpha /= (mesh_.V()*rho_);

    return talpha;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicCloud<CloudType>::rhoEff() const
{
    tmp<volScalarField> trhoEff
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":rhoEff",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimDensity, 0.0)
        )
    );

    scalarField& rhoEff = trhoEff().internalField();
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();

        rhoEff[cellI] += p.nParticle()*p.mass();
    }

    rhoEff /= mesh_.V();

    return trhoEff;
}

template<class CloudType>
inline const Foam::tmp<Foam::volVectorField>
Foam::KinematicCloud<CloudType>::Ul() const
{
    tmp<volVectorField> tUl
    (
        new volVectorField
        (
            IOobject
            (
                this->name() + "Ul",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedVector("zero", dimensionSet(0,1,-1,0,0,0,0), vector(0.0,0.0,0.0))
        )
    );

    tmp<volScalarField> tMass
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "Mass",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(1,0,0,0,0,0,0), 0.0)
        )
    );

    vectorField& Ul = tUl().internalField();
    scalarField& Mass = tMass().internalField();
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();

        Ul[cellI] += p.mass()*p.U()*p.nParticle();
        Mass[cellI] += p.mass()*p.nParticle();
    }

    Ul /= (Mass+VSMALL);

    return tUl;
}

template<class CloudType>
inline const Foam::tmp<Foam::volSymmTensorField>
Foam::KinematicCloud<CloudType>::UlUl() const
{
    tmp<volSymmTensorField> tUlUl
    (
        new volSymmTensorField
        (
            IOobject
            (
                this->name() + "UlUl",
                this->db().time().timeName(),
this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedSymmTensor("zero", dimensionSet(0,2,-2,0,0,0,0),
             symmTensor(0.0,0.0,0.0,0.0,0.0,0.0))
        )
    );

    tmp<volScalarField> tMass
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "Mass",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(1,0,0,0,0,0,0), 0.0)
        )
    );

    symmTensorField& UlUl = tUlUl().internalField();
    scalarField& Mass = tMass().internalField();
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();

        UlUl[cellI] += p.mass()*sqr(p.U())*p.nParticle();
        Mass[cellI] += p.mass()*p.nParticle();
    }

    UlUl /= (Mass+VSMALL);

    return tUlUl;
}


template<class CloudType>
inline const Foam::tmp<Foam::volVectorField>
Foam::KinematicCloud<CloudType>::Uslip() const
{
    tmp<volVectorField> tUslip
    (
        new volVectorField
        (
            IOobject
            (
                this->name() + "Uslip",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedVector("zero", dimensionSet(0,1,-1,0,0,0,0), vector(0.0,0.0,0.0))
        )
    );

tmp<volScalarField> tMass
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "Mass",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(1,0,0,0,0,0,0), 0.0)
        )
    );

    vectorField& Uslip = tUslip().internalField();
    scalarField& Mass = tMass().internalField();
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();
        Uslip[cellI] += p.mass()*(p.Uc()-p.U())*p.nParticle();
        Mass[cellI] += p.mass()*p.nParticle();
    }

    Uslip /= (Mass+VSMALL);

    return tUslip;
}

template<class CloudType>
inline const Foam::tmp<Foam::volSymmTensorField>
Foam::KinematicCloud<CloudType>::UUslip() const                                                                                                           
{
    tmp<volSymmTensorField> tUUslip                                                                                                                       
    (                                                                                                                                                     
        new volSymmTensorField                                                                                                                            
        (                                                                                                                                                 
            IOobject                                                                                                                                      
            (                                                                                                                                             
                this->name() + "UUslip",                                                                                                                  
                this->db().time().timeName(),                                                                                                             
                this->db(),                                                                                                                               
                IOobject::NO_READ,                                                                                                                        
                IOobject::NO_WRITE,                                                                                                                       
                false                                                                                                                                     
            ),                                                                                                                                            
            mesh_,                                                                                                                                        
            dimensionedSymmTensor("zero", dimensionSet(0,2,-2,0,0,0,0),  
symmTensor(0.0,0.0,0.0,0.0,0.0,0.0))                                                                                                         
        )                                                                                                                                                 
    );                                                                                                                                                    


    tmp<volScalarField> tMass                                                                                                                             
    (                                                                                                                                                     
        new volScalarField                                                                                                                                
        (                                                                                                                                                 
            IOobject                                                                                                                                      
            (                                                                                                                                             
                this->name() + "Mass",                                                                                                                    
                this->db().time().timeName(),                                                                                                             
                this->db(),                                                                                                                               
                IOobject::NO_READ,                                                                                                                        
                IOobject::NO_WRITE,                                                                                                                       
                false                                                                                                                                     
            ),                                                                                                                                            
            mesh_,                                                                                                                                        
            dimensionedScalar("zero", dimensionSet(1,0,0,0,0,0,0), 0.0)                                                                                   
        )                                                                                                                                                 
    );                                                                                                                                                    


    symmTensorField& UUslip = tUUslip().internalField();                                                                                                  
    scalarField& Mass = tMass().internalField();                                                                                                          
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)                                                                                      
    {                                                                                                                                                     
        const parcelType& p = iter();                                                                                                                     
        const label cellI = p.cell();                                                                                                                     


        UUslip[cellI] += p.mass()*sqr(p.Uc()-p.U())*p.nParticle();                                                                                        
        Mass[cellI] += p.mass()*p.nParticle();                                                                                                            
    }                                                                                                                                                     


    UUslip /= (Mass+VSMALL); 
return tUUslip;                                                                                                                                       
}

template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicCloud<CloudType>::Udisp() const
{
    tmp<volScalarField> tUdisp
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "Udisp",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(0,0,0,0,0,0,0), 0.0)
        )
    );

    tmp<volScalarField> tMass
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "Mass",
                this->db().time().timeName(),                                                                                                             
                this->db(),                                                                                                                               
                IOobject::NO_READ,                                                                                                                        
                IOobject::NO_WRITE,                                                                                                                       
                false                                                                                                                                     
            ),                                                                                                                                            
            mesh_,                                                                                                                                        
            dimensionedScalar("zero", dimensionSet(1,0,0,0,0,0,0), 0.0)                                                                                   
        )                                                                                                                                                 
    );                                                                                                                                                    
                                                                                                                                                          
    scalarField& Udisp = tUdisp().internalField();                                                                                                        
    scalarField& Mass = tMass().internalField();                                                                                                          
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)                                                                                      
    {                                                                                                                                                     
        const parcelType& p = iter();                                                                                                                     
        const label cellI = p.cell();
        Udisp[cellI] += p.mass()*mag(p.UTurb())/
                     (mag(p.Uc()-p.UTurb()-p.U())+mag(p.UTurb())+VSMALL)*p.nParticle();
Mass[cellI] += p.mass()*p.nParticle();
    }                                                                                                                                                     


    Udisp /= (Mass+VSMALL);                                                                                                                               


    return tUdisp;                                                                                                                                        
}

//Chi-Wei Tsang (05/2012), calculate droplets drag force which is used for
// source term in ksgs equation
template<class CloudType>
inline const Foam::tmp<Foam::volVectorField>
Foam::KinematicCloud<CloudType>::dragForce() const
{

    tmp<volVectorField> tdrag
    (
        new volVectorField
        (
            IOobject
            (
                this->name() + "Drag",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedVector("zero" , dimensionSet(1,-2,-2,0,0,0,0), vector(0.0,0.0,0.0))
        )

    );

    vectorField& drag = tdrag().internalField();
    drag = -UTransDrag().field();

    /*forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();
        const vector Ucc = p.Uc();
        const scalar Urel = mag(p.Uc() - p.U());
        const scalar Re = p.d()*p.rhoc()*Urel/p.muc();
              scalar Cd;

        //eq.(2.39), Nidheesh's Phd thesis
        if (Re <= 1000.0)
        {
            Cd = 24.0/(Re+VSMALL)*(1.0+1.0/6.0*pow((Re+VSMALL), 2.0/3.0));
        }
        else
        {
            Cd = 0.424;
        }

        //eq.(2.38), Nidheesh's PhD thesis
        drag[cellI] += (3.0/8.0*p.rhoc()/(p.rho()+VSMALL)*p.mass()*Urel/(p.d()+VSMALL)*2*
                       (p.Uc() - p.U())*Cd) * p.nParticle();
        //drag[cellI] += p.drag();
    }*/

    drag /= (mesh_.V()*this->db().time().deltaT());
 //drag /= mesh_.V();

    return tdrag;

}

template<class CloudType>
inline const Foam::tmp<Foam::volTensorField>
Foam::KinematicCloud<CloudType>::dipole() const
{
    tmp<volTensorField> tDip
    (
        new volTensorField
        (
            IOobject
            (
                this->name() + "tDip",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedTensor("zero" , dimensionSet(1,-1,-2,0,0,0,0),
            tensor(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0))
       )
   );

   tensorField& Dip = tDip().internalField();
   Dip = DiTrans().field();

   Dip /= (mesh_.V()*this->db().time().deltaT());

   return tDip;
}

//momentum source term due to vaporization, Chi-Wei Tsang (05/2013)
template<class CloudType>
inline const Foam::tmp<Foam::volVectorField>
Foam::KinematicCloud<CloudType>::evap() const
{

    tmp<volVectorField> tevap
    (
        new volVectorField
        (
            IOobject
            (
                this->name() + "Evap",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedVector("zero" , dimensionSet(1,-2,-2,0,0,0,0), vector(0.0,0.0,0.0))
        )

    );

    vectorField& evap = tevap().internalField();
    evap = -(UTrans().field()-UTransDrag().field());

    evap /= (mesh_.V()*this->db().time().deltaT());         
 return tevap;
}

template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicCloud<CloudType>::A_int() const
{

    const scalar pi = 3.1415926;
    tmp<volScalarField> tA_int
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "A_int",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(0,2,0,0,0,0,0), scalar(0.0)),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    volScalarField& A_int = tA_int();
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();

        A_int[cellI] += p.nParticle()*pi*p.d()*p.d();
    }

    //theta.internalField() /= mesh_.V();                                                                                                                 
    A_int.correctBoundaryConditions();                                                                                                                    


    return tA_int;                                                                                                                                        
}

template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicCloud<CloudType>::Wsgs() const
{

    tmp<volScalarField> tWsgs
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "Wsgs",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(1,-1,-3,0,0,0,0), scalar(0.0)),
            zeroGradientFvPatchScalarField::typeName
)
    );

    scalarField& Wsgs = tWsgs().internalField();
    Wsgs = -ksgsTrans().field();

    /*forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();

        Wsgs[cellI] += p.drag() & p.UTurb();
    }*/

    Wsgs /= (mesh_.V()*this->db().time().deltaT());
    //Wsgs.correctBoundaryConditions();



    return tWsgs;

}

// ************************************************************************* //

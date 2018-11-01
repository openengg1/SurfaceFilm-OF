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

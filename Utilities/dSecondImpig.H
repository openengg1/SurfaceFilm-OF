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







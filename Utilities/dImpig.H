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






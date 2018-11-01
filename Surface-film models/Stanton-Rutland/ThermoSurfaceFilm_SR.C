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

#include "ThermoSurfaceFilm.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "Pstream.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class CloudType>
Foam::wordList Foam::ThermoSurfaceFilm<CloudType>::interactionTypeNames_
(
    IStringStream
    (
        "(absorb bounce splashBai)"
    )()
);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class CloudType>
typename Foam::ThermoSurfaceFilm<CloudType>::interactionType
Foam::ThermoSurfaceFilm<CloudType>::interactionTypeEnum(const word& it) const
{
    forAll(interactionTypeNames_, i)
    {
        if (interactionTypeNames_[i] == it)
        {
            return interactionType(i);
        }
    }

    FatalErrorIn
    (
        "ThermoSurfaceFilm<CloudType>::interactionType "
        "ThermoSurfaceFilm<CloudType>::interactionTypeEnum"
        "("
            "const word& it"
        ") const"
    )   << "Unknown interaction type " << it
        << ". Valid interaction types include: " << interactionTypeNames_
        << abort(FatalError);

    return interactionType(0);
}


template<class CloudType>
Foam::word Foam::ThermoSurfaceFilm<CloudType>::interactionTypeStr
(
    const interactionType& it
) const
{
    if (it >= interactionTypeNames_.size())
    {
        FatalErrorIn
        (
            "ThermoSurfaceFilm<CloudType>::interactionType "
            "ThermoSurfaceFilm<CloudType>::interactionTypeStr"
            "("
                "const interactionType& it"
            ") const"
        )   << "Unknown interaction type enumeration" << abort(FatalError);
    }

    return interactionTypeNames_[it];
}


template<class CloudType>
Foam::vector Foam::ThermoSurfaceFilm<CloudType>::tangentVector
(
    const vector& v
) const
{
    vector tangent = vector::zero;
    scalar magTangent = 0.0;

    while (magTangent < SMALL)
    {
        vector vTest = rndGen_.sample01<vector>();
        tangent = vTest - (vTest & v)*v;
        magTangent = mag(tangent);
    }

    return tangent/magTangent;
}

//Rohit Mishra (07/10/2018)
//Stanton-Rutland Weibull distribution
template<class CloudType>
Foam::vector Foam::ThermoSurfaceFilm<CloudType>::splashDirection
(
    const vector Ut,
    const vector Un,
    const scalar theta
) const
{
scalar minValue = 0.0;
scalar maxValue = 1.0;
scalar d = 0.158*exp(0.017*theta);
scalar n = 1.1 + 0.02*theta;

if (theta <= 50)
{
n=2.1;
}

else
{
//Default condition
}
scalar K = 1.0 - exp(-pow((maxValue - minValue)/d, n));
scalar y = rndGen_.sample01<scalar>();
scalar x = minValue + d*::pow(-log(1.0 - y*K), 1.0/n);
vector Vn=-x*Un;
scalar theta_s=65.4+0.266*(90-theta);  
theta_s=pi/180*(90-theta_s);
scalar magVt=mag(Vn)/tan(theta_s);

//Rohit Mishra (07/16/2018)
//Naber-Reitz Azimuthal equation
scalar alpha=90-theta;
scalar Err=10;
scalar beta=10;
scalar diff=1;
scalar alpha_p=90;
while (Err>1e-3)
{ 
   y=((exp(beta)+1)/(exp(beta)-1))*(1/(1+pow((pi/beta),2)));
    alpha_p=asin(y )*180/pi;
    Err=abs(alpha_p-alpha);
    if (Err<40)
    {
        if (alpha>=0 && alpha <= 1)
        { 
         diff=Err/100;
        }
        if (alpha>1 && alpha<=70)
        {            
         diff=Err/10;
        }
        if (alpha>70 && alpha<=90)
        {
          diff=Err;
        }
    }
    if (alpha_p<alpha)
    { 
     Err=1e-4;
    }
    
    beta=beta-diff;
}

beta=beta+diff;
scalar p=rndGen_.sample01<scalar>();
scalar psi=(-pi/beta) * log (1 - p * (1-exp(-beta)));

vector xvec(1,0,0);
scalar gamma=acos((Ut & xvec)/(mag(Ut)*mag(xvec)));
scalar sign=1;

if (p<0.5)
{
sign=-1;
}

scalar Vtz=1;
scalar Vtx=Vtz/tan(gamma-sign*psi);
vector Vt(Vtx,0,Vtz);
Vt=Vt/mag(Vt);
Vt=Vt*magVt;

Info<<endl<<"Ut="<<Ut<<endl<<"Un="<<Un<<endl;
Info<<endl<<"gamma="<<gamma*180/pi<<endl;
Info<<endl<<"Beta="<<beta<<endl;
Info<<endl<<"Azimuthal angle="<<psi*180/pi<<endl;

Info<<endl<<"Impingement angle="<<theta<<endl<<"Reflection angle="<<(180/pi)*(theta_s)<<endl<<"Vt="<<Vt<<endl<<"Vn="<<Vn<<endl<<"Weibull X="<<x<<endl;
return (Vt+Vn)/mag(Vt+Vn);

//Info<<endl<<"Impingement angle="<<theta<<endl<<"Reflection angle="<<theta_s<<endl;


/*
   // azimuthal angle [rad]
    const scalar phiSi = twoPi*rndGen_.sample01<scalar>();

// ejection angle [rad]
    scalar thetaSi = pi/180.0*(rndGen_.sample01<scalar>()*(50 - 5));

// Rohit Mishra (07/03/2018)


if (sp_U==nf)
{
    // ejection angle [rad]
Info<<endl<<endl<<"Normal impact"<<endl<<endl;
    thetaSi = pi/180.0*(rndGen_.sample01<scalar>()*(50 - 5) + 5);
}
else
{
Info<<endl<<endl<<"Oblique impact"<<endl<<endl;
// Default condition
}


    // direction vector of new parcel
    const scalar alpha = sin(thetaSi);
    const scalar dcorr = cos(thetaSi);
    const vector normal = alpha*(tanVec1*cos(phiSi) + tanVec2*sin(phiSi));
    vector dirVec = dcorr*sp_U;
    dirVec += normal;

    return dirVec/mag(dirVec);
*/



}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::absorbInteraction
(
    regionModels::surfaceFilmModels::surfaceFilmModel& filmModel,
    const parcelType& p,
    const polyPatch& pp,
    const label faceI,
    const scalar mass,
    bool& keepParticle
)
{
    if (debug)
    {
        Info<< "Parcel " << p.origId() << " absorbInteraction" << endl;
    }

    // Patch face normal
    const vector& nf = pp.faceNormals()[faceI];

    // Patch velocity
    const vector& Up = this->owner().U().boundaryField()[pp.index()][faceI];

    // Relative parcel velocity
    const vector Urel = p.U() - Up;

    // Parcel normal velocity
    const vector Un = nf*(Urel & nf);

    // Parcel tangential velocity
    const vector Ut = Urel - Un;

    filmModel.addSources
    (
        pp.index(),
        faceI,
        mass,                           // mass
        mass*Ut,                        // tangential momentum
        mass*mag(Un),                   // impingement pressure
        mass*p.hs()                     // energy
    );

    this->nParcelsTransferred()++;

    keepParticle = false;
}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::bounceInteraction
(
    parcelType& p,
    const polyPatch& pp,
    const label faceI,
    bool& keepParticle
) const
{
    if (debug)
    {
        Info<< "Parcel " << p.origId() << " bounceInteraction" << endl;
    }

    // Patch face normal
    const vector& nf = pp.faceNormals()[faceI];

    // Patch velocity
    const vector& Up = this->owner().U().boundaryField()[pp.index()][faceI];

    // Relative parcel velocity
    const vector Urel = p.U() - Up;

    // Flip parcel normal velocity component
    p.U() -= 2.0*nf*(Urel & nf);

    keepParticle = true;
}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::drySplashInteraction
(
    regionModels::surfaceFilmModels::surfaceFilmModel& filmModel,
    const parcelType& p,
    const polyPatch& pp,
    const label faceI,
    bool& keepParticle
)
{
    if (debug)
    {
        Info<< "Parcel " << p.origId() << " drySplashInteraction" << endl;
    }

    const liquidProperties& liq = thermo_.liquids().properties()[0];

    // Patch face velocity and normal
    const vector& Up = this->owner().U().boundaryField()[pp.index()][faceI];
    const vector& nf = pp.faceNormals()[faceI];

    // local pressure
    const scalar pc = thermo_.thermo().p()[p.cell()];

    // Retrieve parcel properties
    const scalar m = p.mass()*p.nParticle();
    const scalar rho = p.rho();
    const scalar d = p.d();
    const scalar sigma = liq.sigma(pc, p.T());
    const scalar mu = liq.mu(pc, p.T());
    const vector Urel = p.U() - Up;
    const vector Un = nf*(Urel & nf);

    // Laplace number
    const scalar La = rho*sigma*d/sqr(mu);

    // Weber number
    const scalar We = rho*magSqr(Un)*d/sigma;
Info<<endl<<endl<<"rho="<<rho<<endl<<"Diameter of parcel="<<d<<endl<<"Laplacian number="<<La<<endl<<"Viscosity="<<mu<<endl<<"Surface Tension="<<sigma<<endl;

    // Critical Weber number
    const scalar Wec = Adry_*pow(La, -0.183);
Info<<endl<<"We="<<We<<endl;
Info<<endl<<"Wec="<<Wec<<endl;

    if (We < Wec) // adhesion - assume absorb
    {
        absorbInteraction(filmModel, p, pp, faceI, m, keepParticle);
    }
    else // splash
    {
       // ratio of incident mass to splashing mass
    //   const scalar mRatio = 0.2 + 0.6*rndGen_.sample01<scalar>();




//Rohit Mishra (08/24/2018)
//Stanton-Rutland mass ratio curve fit 1996

scalar f = mag(Un)/d;
scalar Udim=mag(Un)*(pow((rho/sigma),(1/4)))*pow(mu,(-1/8))*pow(f,(-3/8));
scalar mRatio=-27.2 + 3.15*mag(Udim) - 0.116*pow(mag(Udim),2) + 1.4*1e-3*pow(mag(Udim),3);

label nTrans0 =
        this->template getModelProperty<label>("nParcelsTransferred");

label nTransTotal =
        nTrans0 + this->nParcelsTransferred();

label nSplash0 = this->template getModelProperty<label>("nParcelsSplashed");

label nSplashTotal =
        nSplash0 + returnReduce(nParcelsSplashed_, sumOp<label>());



if ((nTransTotal+nSplashTotal)>1)
{
f=(nTransTotal+nSplashTotal)/this->owner().db().time().value();

Udim=mag(Un)*(pow((rho/sigma),(1/4)))*pow(mu,(-1/8))*pow(f,(-3/8));


mRatio = -27.2 + 3.15*mag(Udim) - 0.116*pow(mag(Udim),2) + 1.4*1e-3*pow(mag(Udim),3);

}

if (mRatio > 1)
{
mRatio = 0.71;
}

Info<<endl<<"Udim:"<<Udim<<endl;
Info<<endl<<"Normal velocity="<<mag(Un)<<endl;
Info<<endl<<"nTransTotal="<<nTransTotal;
Info<<endl<<"nSplashTotal="<<nSplashTotal;
Info<<endl<<"Current time="<<this->owner().db().time().value()<<endl;
Info<<endl<<"Frequency="<<f<<endl<<"Mass Ratio="<<mRatio<<endl;

        splashInteraction
            (filmModel, p, pp, faceI, mRatio, We, Wec, sigma, keepParticle);
  


  }
}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::wetSplashInteraction
(
    regionModels::surfaceFilmModels::surfaceFilmModel& filmModel,
    parcelType& p,
    const polyPatch& pp,
    const label faceI,
    bool& keepParticle
)
{
    if (debug)
    {
        Info<< "Parcel " << p.origId() << " wetSplashInteraction" << endl;
    }

    const liquidProperties& liq = thermo_.liquids().properties()[0];

    // Patch face velocity and normal
    const vector& Up = this->owner().U().boundaryField()[pp.index()][faceI];
    const vector& nf = pp.faceNormals()[faceI];

    // local pressure
    const scalar pc = thermo_.thermo().p()[p.cell()];

    // Retrieve parcel properties
    const scalar m = p.mass()*p.nParticle();
    const scalar rho = p.rho();
    const scalar d = p.d();
    vector& U = p.U();
    const scalar sigma = liq.sigma(pc, p.T());
    const scalar mu = liq.mu(pc, p.T());
    const vector Urel = p.U() - Up;
    const vector Un = nf*(Urel & nf);
    const vector Ut = Urel - Un;

    // Laplace number
    const scalar La = rho*sigma*d/sqr(mu);

    // Weber number
    const scalar We = rho*magSqr(Un)*d/sigma;

    // Critical Weber number
  // const scalar Wec = Awet_*pow(La, -0.183);
//Info << endl<<"Normal velocity="<<mag(Un)<<endl;
Info<<endl<<endl<<"rho="<<rho<<endl<<"Diameter of parcel="<<d<<endl<<"Laplacian number="<<La<<endl<<"Viscosity="<<mu<<endl<<"Surface Tension="<<sigma<<endl;

   

scalar f = mag(Un)/35e-3;

label nTrans0 =
        this->template getModelProperty<label>("nParcelsTransferred");

label nTransTotal =
        nTrans0 + this->nParcelsTransferred();

label nSplash0 = this->template getModelProperty<label>("nParcelsSplashed");

label nSplashTotal =
        nSplash0 + returnReduce(nParcelsSplashed_, sumOp<label>());



if ((nTransTotal+nSplashTotal)>1)
{
f=(nTransTotal+nSplashTotal)/this->owner().db().time().value();
}


Info<<endl << "Parameter:"<<"Dia:"<<d<<"rho:"<<rho<<"sigma:"<<sigma<<"mu:"<<mu<<"f:"<<f<<endl;

Info << endl << "Power:"<<"pow(18,2):"<<pow(18,2)<<" pow((rho/sigma),(1/2)):"<<pow((rho/sigma),0.5)<<" pow(mu,(1/4)):"<<pow(mu,0.25)<<" pow(f,0.75):"<<pow(f,0.75)<<endl;
const scalar Wec = pow(18,2) * d * pow((rho/sigma),0.5) * pow(mu,0.25) * pow(f,0.75);

Info<< endl << " We = "<<We<<endl<<" Wec = " <<Wec<<endl;

    if (We < 5) // adhesion - assume absorb
    {
        absorbInteraction(filmModel, p, pp, faceI, m, keepParticle);
    }
    else if ((We >= 5) && (We < 10)) // bounce
    {
        // incident angle of impingement
        const scalar theta = pi/2 - acos(U/mag(U) & nf);

        // restitution coefficient
        const scalar epsilon = 0.993 - theta*(1.76 - theta*(1.56 - theta*0.49));

        // update parcel velocity
        U = -epsilon*(Un) + 5/7*(Ut);

        keepParticle = true;
        return;
    }
    else if ((We >= 10) && (We < Wec)) // spread - assume absorb
    {
        absorbInteraction(filmModel, p, pp, faceI, m, keepParticle);
    }
    else    // splash
    {
        // ratio of incident mass to splashing mass
       // splash mass can be > incident mass due to film entrainment
       // const scalar mRatio = 0.2 + 0.9*rndGen_.sample01<scalar>();

//Rohit Mishra (08/24/2018)
//Stanton-Rutland mass ratio curve fit 1996

scalar Udim=mag(Un)*(pow((rho/sigma),(1/4)))*pow(mu,(-1/8))*pow(f,(-3/8));
scalar mRatio=-27.2 + 3.15*mag(Udim) - 0.116*pow(mag(Udim),2) + 1.4*1e-3*pow(mag(Udim),3);

if(mRatio < 0)
{
mRatio=0.001;
}

if (mRatio > 1)
{
mRatio = 0.71;
}


Info<<endl<<"Udim:"<<Udim<<endl;
Info<<endl<<"Normal velocity="<<mag(Un)<<endl;
Info<<endl<<"nTransTotal="<<nTransTotal;
Info<<endl<<"nSplashTotal="<<nSplashTotal;
Info<<endl<<"Current time="<<this->owner().db().time().value()<<endl;
Info<<endl<<"Frequency="<<f<<endl<<"Mass Ratio="<<mRatio<<endl;

  
     splashInteraction
            (filmModel, p, pp, faceI, mRatio, We, Wec, sigma, keepParticle,f);
    }
}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::splashInteraction
(
    regionModels::surfaceFilmModels::surfaceFilmModel& filmModel,
    const parcelType& p,
    const polyPatch& pp,
    const label faceI,
    const scalar mRatio,
    const scalar We,
    const scalar Wec,
    const scalar sigma,
    bool& keepParticle,
    scalar f
)
{


 const liquidProperties& liq = thermo_.liquids().properties()[0];
  const scalar pc = thermo_.thermo().p()[p.cell()];
 
  // Patch face velocity and normal
    const fvMesh& mesh = this->owner().mesh();
    const vector& Up = this->owner().U().boundaryField()[pp.index()][faceI];
    const vector& nf = pp.faceNormals()[faceI];

    // Determine direction vectors tangential to patch normal
    const vector tanVec1 = tangentVector(nf);
    const vector tanVec2 = nf^tanVec1;

    // Retrieve parcel properties
    const scalar np = p.nParticle();
    const scalar m = p.mass()*np;
    const scalar d = p.d();
    const vector Urel = p.U() - Up;
    const vector Un = nf*(Urel & nf);
    const vector Ut = Urel - Un;
    const scalar rho = p.rho();
    const scalar mu = liq.mu(pc, p.T());
    


    const vector& posC = mesh.C()[p.cell()];
    const vector& posCf = mesh.Cf().boundaryField()[pp.index()][faceI];
   //vector& U = p.U();    
       // incident angle of impingement
        const scalar theta = 180/pi*(pi/2 - acos(p.U()/mag(p.U()) & nf));



    // total mass of (all) splashed parcels
    const scalar mSplash = m*mRatio;

Info<<endl<<"Incoming mass=" <<m<<endl<<"Splashed mass="<<mSplash<<endl;
Info<<endl<<"Incoming velocity="<<mag(Urel)<<endl;
Info<< endl << " We = "<<We<<endl<<" Wec = " <<Wec<<endl;

scalar b = 2.71 - 9.25e-4*We;
scalar theta_w = 0.21 - 7.69e-5*We;

    // cumulative diameter splash distribution
    const scalar dMax = 0.9*cbrt(mRatio)*d;
    const scalar dMin = 0.1*dMax;
    const scalar K =  1.0 - exp(-pow((dMax - dMin)/(b*d), theta_w));

    // surface energy of secondary parcels [J]
    scalar ESigmaSec = 0;

    // sample splash distribution to determine secondary parcel diameters
    scalarList dNew(parcelsPerSplash_);
    scalarList npNew(parcelsPerSplash_);


Info<<endl<<"Incoming diamter="<<d<<endl;
    forAll(dNew, i)
    {
        scalar y = rndGen_.sample01<scalar>();
        scalar x = (dMin/d) + b*::pow(-log(1.0 - y*K), 1.0/theta_w);
        dNew[i] = x*d;
        npNew[i] = mRatio*np*pow3(d)/pow3(dNew[i])/parcelsPerSplash_;
        ESigmaSec += npNew[i]*sigma*p.areaS(dNew[i]);
Info << endl << "Splashed diameter="<<dNew[i]<<endl<< "Splashed particles="<<npNew[i]<<endl;



    }



    // incident kinetic energy [J]
    const scalar EKIn = 0.5*m*magSqr(Urel);

    // incident surface energy [J]
    const scalar ESigmaIn = np*sigma*p.areaS(d);
    
    //Velocity at splashing threshold
    const scalar Vs = 18 * pow((sigma/rho),1/4) * pow(mu,1/8) * pow(f,3/8);

    // dissipative energy
    const scalar Ed = 0.5 * m * sqr(Vs);

    // total energy [J]
    const scalar EKs = EKIn + ESigmaIn - ESigmaSec - Ed;

    // switch to absorb if insufficient energy for splash
    if (EKs <= 0)
    {
        absorbInteraction(filmModel, p, pp, faceI, m, keepParticle);
        return;
    }

    // helper variables to calculate magUns0
    const scalar logD = log(d);
    const scalar coeff2 = log(dNew[0]) - logD + ROOTVSMALL;
    scalar coeff1 = 0.0;
    forAll(dNew, i)
    {
        coeff1 += sqr(log(dNew[i]) - logD);
    }

    // magnitude of the normal velocity of the first splashed parcel
    const scalar magUns0 =
        sqrt(2.0*parcelsPerSplash_*EKs/mSplash/(1.0 + coeff1/sqr(coeff2)));

    // Set splashed parcel properties
    forAll(dNew, i)
    {
        const vector dirVec = splashDirection(Ut,Un,theta);

        // Create a new parcel by copying source parcel
        parcelType* pPtr = new parcelType(p);

        pPtr->origId() = pPtr->getNewParticleID();

        pPtr->origProc() = Pstream::myProcNo();

        if (splashParcelType_ >= 0)
        {
            pPtr->typeId() = splashParcelType_;
        }

        // perturb new parcels towards the owner cell centre
        pPtr->position() += 0.5*rndGen_.sample01<scalar>()*(posC - posCf);

        pPtr->nParticle() = npNew[i];

        pPtr->d() = dNew[i];

        pPtr->U() = dirVec*(mag(Cf_*Ut) + magUns0*(log(dNew[i]) - logD)/coeff2);

        Info << endl << "Splashed velocity="<< (mag(Cf_*Ut) + magUns0*(log(dNew[i]) - logD)/coeff2)<<endl;

        // Apply correction to velocity for 2-D cases
        meshTools::constrainDirection(mesh, mesh.solutionD(), pPtr->U());

        // Add the new parcel
        this->owner().addParticle(pPtr);

        nParcelsSplashed_++;
    }

    // transfer remaining part of parcel to film 0 - splashMass can be -ve
    // if entraining from the film
    const scalar mDash = m - mSplash;
    absorbInteraction(filmModel, p, pp, faceI, mDash, keepParticle);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermoSurfaceFilm<CloudType>::ThermoSurfaceFilm
(
    const dictionary& dict,
    CloudType& owner
)
:
    SurfaceFilmModel<CloudType>(dict, owner, typeName),
    rndGen_(owner.rndGen()),
    thermo_
    (
        owner.db().objectRegistry::template lookupObject<SLGThermo>("SLGThermo")
    ),
    TFilmPatch_(0),
    CpFilmPatch_(0),
    interactionType_
    (
        interactionTypeEnum(this->coeffDict().lookup("interactionType"))
    ),
    h0_(0.0),
    deltaWet_(0.0),
    splashParcelType_(0),
    parcelsPerSplash_(0),
    Adry_(0.0),
    Awet_(0.0),
    Cf_(0.0),
    nParcelsSplashed_(0)
{
    Info<< "    Applying " << interactionTypeStr(interactionType_)
        << " interaction model" << endl;

    if (interactionType_ == itSplashBai)
    {
        this->coeffDict().lookup("deltaWet") >> deltaWet_;
        this->coeffDict().lookup("h0") >> h0_;
        splashParcelType_ =
            this->coeffDict().lookupOrDefault("splashParcelType", -1);
        parcelsPerSplash_ =
            this->coeffDict().lookupOrDefault("parcelsPerSplash", 2);
        this->coeffDict().lookup("Adry") >> Adry_;
        this->coeffDict().lookup("Awet") >> Awet_;
        this->coeffDict().lookup("Cf") >> Cf_;
    }
}


template<class CloudType>
Foam::ThermoSurfaceFilm<CloudType>::ThermoSurfaceFilm
(
    const ThermoSurfaceFilm<CloudType>& sfm
)
:
    SurfaceFilmModel<CloudType>(sfm),
    rndGen_(sfm.rndGen_),
    thermo_(sfm.thermo_),
    TFilmPatch_(sfm.TFilmPatch_),
    CpFilmPatch_(sfm.CpFilmPatch_),
    interactionType_(sfm.interactionType_),
    h0_(sfm.h0_),
    deltaWet_(sfm.deltaWet_),
    splashParcelType_(sfm.splashParcelType_),
    parcelsPerSplash_(sfm.parcelsPerSplash_),
    Adry_(sfm.Adry_),
    Awet_(sfm.Awet_),
    Cf_(sfm.Cf_),
    nParcelsSplashed_(sfm.nParcelsSplashed_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermoSurfaceFilm<CloudType>::~ThermoSurfaceFilm()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::ThermoSurfaceFilm<CloudType>::transferParcel
(
    parcelType& p,
    const polyPatch& pp,
    bool& keepParticle
)
{
    // Retrieve the film model from the owner database
    regionModels::surfaceFilmModels::surfaceFilmModel& filmModel =
        const_cast<regionModels::surfaceFilmModels::surfaceFilmModel&>
        (
            this->owner().db().time().objectRegistry::template
                lookupObject<regionModels::surfaceFilmModels::surfaceFilmModel>
                (
                    "surfaceFilmProperties"
                )
        );

    const label patchI = pp.index();

    if (filmModel.isRegionPatch(patchI))
    {
        const label faceI = pp.whichFace(p.face());

        switch (interactionType_)
        {
            case itBounce:
            {
                bounceInteraction(p, pp, faceI, keepParticle);
Info<<endl<<endl<<"Bounce interaction selected" <<endl<<endl<<endl;
                break;
            }
            case itAbsorb:
            {
                const scalar m = p.nParticle()*p.mass();
                absorbInteraction(filmModel, p, pp, faceI, m, keepParticle);
Info<<endl<<endl<<"absorb interaction selected" <<endl<<endl<<endl;
                break;
            }
            case itSplashBai:
            {
               
Info<<endl<<endl<<"Wet splash interaction selected" <<endl<<endl<<endl;
         wetSplashInteraction(filmModel, p, pp, faceI, keepParticle);
                

                break;
            }
            default:
            {
                FatalErrorIn
                (
                    "bool ThermoSurfaceFilm<CloudType>::transferParcel"
                    "("
                        "parcelType&, "
                        "const polyPatch&, "
                        "bool&"
                    ")"
                )   << "Unknown interaction type enumeration"
                    << abort(FatalError);
            }
        }

        // transfer parcel/parcel interactions complete
        return true;
    }

    // parcel not interacting with film
    return false;
}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::cacheFilmFields
(
    const label filmPatchI,
    const label primaryPatchI,
    const regionModels::surfaceFilmModels::surfaceFilmModel& filmModel
)
{
    SurfaceFilmModel<CloudType>::cacheFilmFields
    (
        filmPatchI,
        primaryPatchI,
        filmModel
    );

    TFilmPatch_ = filmModel.Ts().boundaryField()[filmPatchI];
    filmModel.toPrimary(filmPatchI, TFilmPatch_);

    CpFilmPatch_ = filmModel.Cp().boundaryField()[filmPatchI];
    filmModel.toPrimary(filmPatchI, CpFilmPatch_);
}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::setParcelProperties
(
    parcelType& p,
    const label filmFaceI
) const
{
    SurfaceFilmModel<CloudType>::setParcelProperties(p, filmFaceI);

    // Set parcel properties
    p.T() = TFilmPatch_[filmFaceI];
    p.Cp() = CpFilmPatch_[filmFaceI];
}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::info(Ostream& os)
{
    SurfaceFilmModel<CloudType>::info(os);

    label nSplash0 = this->template getModelProperty<label>("nParcelsSplashed");
    label nSplashTotal =
        nSplash0 + returnReduce(nParcelsSplashed_, sumOp<label>());

    os  << "    New film splash parcels         = " << nSplashTotal << endl;

    if (this->outputTime())
    {
        this->setModelProperty("nParcelsSplashed", nSplashTotal);
        nParcelsSplashed_ = 0;
    }
}


// ************************************************************************* //

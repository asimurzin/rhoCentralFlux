#!/usr/bin/env python

#--------------------------------------------------------------------------------------
## pythonFlu - Python wrapping for OpenFOAM C++ API
## Copyright (C) 2010- Alexey Petrov
## Copyright (C) 2009-2010 Pebble Bed Modular Reactor (Pty) Limited (PBMR)
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 
## See http://sourceforge.net/projects/pythonflu
##
## Author : Alexey PETROV, Andrey Simurzin
##


#---------------------------------------------------------------------------
from BCs import rho
from Foam import ref, man

#---------------------------------------------------------------------------
def _rhoBoundaryTypes( p ):
    pbf = p.ext_boundaryField()
    rhoBoundaryTypes = pbf.types()
    
    for patchi in range( rhoBoundaryTypes.size() ):
        if str( rhoBoundaryTypes[patchi] ) == "waveTransmissive":
           rhoBoundaryTypes[patchi] = ref.zeroGradientFvPatchScalarField.typeName
           pass
        elif pbf[patchi].fixesValue():
           from BCs.rho import fixedRhoFvPatchScalarField
           rhoBoundaryTypes[patchi] = fixedRhoFvPatchScalarField.typeName
           pass
        pass
    
    return pbf, rhoBoundaryTypes


#-------------------------------------------------------------------------------------
def readThermophysicalProperties( runTime, mesh ):
    ref.ext_Info() << "Reading thermophysicalProperties\n" << ref.nl
    
    # Pr defined as a separate constant to enable calculation of k, currently
    # inaccessible through thermo
    thermophysicalProperties = man.IOdictionary( man.IOobject( ref.word( "thermophysicalProperties" ),
                                                               ref.fileName( runTime.constant() ),
                                                               mesh,
                                                               ref.IOobject.MUST_READ,
                                                               ref.IOobject.NO_WRITE ) )
    
    Pr = ref.dimensionedScalar.lookupOrDefault( ref.word( "Pr" ), thermophysicalProperties(), 1.0 )
    
    return thermophysicalProperties, Pr


#---------------------------------------------------------------------------
def _createFields( runTime, mesh ):
    ref.ext_Info() << "Reading thermophysical properties\n" << ref.nl
    
    thermo = man.basicPsiThermo.New( mesh )
    
    p = man.volScalarField( thermo.p(), man.Deps( thermo ) )
    e = man.volScalarField( thermo.e(), man.Deps( thermo ) ) 
    T = man.volScalarField( thermo.T(), man.Deps( thermo ) )
    psi = man.volScalarField( thermo.psi(), man.Deps( thermo ) )
    mu = man.volScalarField( thermo.mu(), man.Deps( thermo ) )
    
    inviscid = True
    if mu.internalField().max() > 0.0:
       inviscid = False
       pass
    
    ref.ext_Info() << "Reading field U\n" << ref.nl
    U = man.volVectorField( man.IOobject( ref.word( "U" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.MUST_READ,
                                          ref.IOobject.AUTO_WRITE ),
                            mesh )
    
    pbf, rhoBoundaryTypes = _rhoBoundaryTypes( p )
    
    rho = man.volScalarField( man.IOobject( ref.word( "rho" ),
                                            ref.fileName( runTime.timeName() ),
                                            mesh,
                                            ref.IOobject.NO_READ,
                                            ref.IOobject.AUTO_WRITE ),
                             man.volScalarField( thermo.rho(), man.Deps( thermo ) ),
                             rhoBoundaryTypes )
    rhoU = man.volVectorField( man.IOobject( ref.word( "rhoU" ),
                                             ref.fileName( runTime.timeName() ),
                                             mesh,
                                             ref.IOobject.NO_READ,
                                             ref.IOobject.NO_WRITE ),
                               rho*U )
    rhoE = man.volScalarField( man.IOobject( ref.word( "rhoE" ),
                                             ref.fileName( runTime.timeName() ),
                                             mesh,
                                             ref.IOobject.NO_READ,
                                             ref.IOobject.NO_WRITE ),
                               rho * ( e + man.volScalarField( 0.5 * U.magSqr(), man.Deps( U ) ) ) )
    
    pos = man.surfaceScalarField( man.IOobject( ref.word( "pos" ),
                                                ref.fileName( runTime.timeName() ),
                                                mesh ),
                                  mesh,
                                  ref.dimensionedScalar( ref.word( "pos" ), ref.dimless, 1.0) )
    
    neg = man.surfaceScalarField( man.IOobject( ref.word( "neg" ),
                                                ref.fileName( runTime.timeName() ),
                                                mesh ),
                                  mesh,
                                  ref.dimensionedScalar( ref.word( "neg" ), ref.dimless, -1.0 ) )

   
    phi = man.surfaceScalarField( ref.word( "phi" ),
                                  man.surfaceVectorField( mesh.Sf(), man.Deps( mesh ) ) & man.fvc.interpolate( rhoU ) )
  
    ref.ext_Info() << "Creating turbulence model\n" << ref.nl
    turbulence = man.compressible.turbulenceModel.New( rho, U, phi, thermo )
    
    return thermo, p, e, T, psi, mu, U, pbf, rhoBoundaryTypes, rho, rhoU, rhoE, pos, neg, inviscid, phi, turbulence


#--------------------------------------------------------------------------------------
def readFluxScheme( mesh ):
    fluxScheme = ref.word("Kurganov")
    
    if mesh.schemesDict().readIfPresent( ref.word( "fluxScheme" ), fluxScheme ):
       if str( fluxScheme ) == "Tadmor" or str( fluxScheme ) == "Kurganov":
          ref.ext_Info() << "fluxScheme: " << fluxScheme << ref.nl
          pass
       else:
          ref.ext_Info() << "rhoCentralFoam::readFluxScheme" \
                     << "fluxScheme: " << fluxScheme \
                     << " is not a valid choice. " \
                     << "Options are: Tadmor, Kurganov" <<ref.nl
          import os
          os.abort()           
          pass
       pass
    return fluxScheme

#--------------------------------------------------------------------------------------
def compressibleCourantNo( mesh, amaxSf, runTime ):
    CoNum = 0.0
    meanCoNum = 0.0
    
    if mesh.nInternalFaces():
       amaxSfbyDelta = mesh.deltaCoeffs() * amaxSf
       CoNum = ( amaxSfbyDelta / mesh.magSf() ).ext_max().value() * runTime.deltaT().value()
       meanCoNum = ( amaxSfbyDelta.sum() / mesh.magSf().sum() ).value() * runTime.deltaTValue()
       pass
    ref.ext_Info() << "Mean and max Courant Numbers = " << meanCoNum << " " << CoNum << ref.nl
    
    return CoNum, meanCoNum


#--------------------------------------------------------------------------------------
def main_standalone( argc, argv ):

    args = ref.setRootCase( argc, argv )

    runTime = man.createTime( args )

    mesh = man.createMesh( runTime )
    
    thermo, p, e, T, psi, mu, U, pbf, rhoBoundaryTypes, rho, rhoU, rhoE, pos, \
                          neg, inviscid, phi, turbulence = _createFields( runTime, mesh )
    
    thermophysicalProperties, Pr = readThermophysicalProperties( runTime, mesh )
    
    fluxScheme = readFluxScheme( mesh )
    
    v_zero = ref.dimensionedScalar( ref.word( "v_zero" ), ref.dimVolume / ref.dimTime, 0.0)
    
    ref.ext_Info() << "\nStarting time loopqqqqq\n" << ref.nl
    
    while runTime.run() :
        # --- upwind interpolation of primitive fields on faces
        rho_pos = ref.fvc.interpolate( rho, pos, ref.word( "reconstruct(rho)" ) )
        rho_neg = ref.fvc.interpolate( rho, neg, ref.word( "reconstruct(rho)" ) )
        
        rhoU_pos = ref.fvc.interpolate( rhoU, pos, ref.word( "reconstruct(U)" ) )
        rhoU_neg = ref.fvc.interpolate( rhoU, neg, ref.word( "reconstruct(U)" ) )

        rPsi = 1.0 / psi
        rPsi_pos = ref.fvc.interpolate( rPsi, pos, ref.word( "reconstruct(T)" ) )
        rPsi_neg = ref.fvc.interpolate( rPsi, neg, ref.word( "reconstruct(T)" ) )

        e_pos = ref.fvc.interpolate( e, pos, ref.word( "reconstruct(T)" ) )
        e_neg = ref.fvc.interpolate( e, neg, ref.word( "reconstruct(T)" ) )

        U_pos = rhoU_pos / rho_pos
        U_neg = rhoU_neg / rho_neg

        p_pos = rho_pos * rPsi_pos
        p_neg = rho_neg * rPsi_neg

        phiv_pos = U_pos & mesh.Sf()
        phiv_neg = U_neg & mesh.Sf()

        c = ( thermo.Cp() / thermo.Cv() * rPsi ).sqrt()
        cSf_pos = ref.fvc.interpolate( c, pos, ref.word( "reconstruct(T)" ) ) * mesh.magSf()
        cSf_neg = ref.fvc.interpolate( c, neg, ref.word( "reconstruct(T)" ) ) * mesh.magSf()
   
        ap = ( phiv_pos + cSf_pos ).ext_max( phiv_neg + cSf_neg ).ext_max( v_zero )
        am = ( phiv_pos - cSf_pos ).ext_min( phiv_neg - cSf_neg ).ext_min( v_zero )

        a_pos = ap / ( ap - am )
        
        amaxSf = ref.surfaceScalarField( ref.word( "amaxSf" ), am.mag().ext_max( ap.mag() ) )
        
        aSf = am * a_pos

        if str( fluxScheme ) == "Tadmor":
           aSf << -0.5 * amaxSf
           a_pos << 0.5
           pass

        a_neg = 1.0 - a_pos
        
        phiv_pos *= a_pos
        phiv_neg *= a_neg
        
        aphiv_pos = phiv_pos - aSf
        aphiv_neg = phiv_neg + aSf
        
        # Reuse amaxSf for the maximum positive and negative fluxes
        # estimated by the central scheme
        amaxSf << aphiv_pos.mag().ext_max(  aphiv_neg.mag() )

        CoNum, meanCoNum = compressibleCourantNo( mesh, amaxSf, runTime )
        
        adjustTimeStep, maxCo, maxDeltaT = ref.readTimeControls( runTime )
        
        runTime = ref.setDeltaT( runTime, adjustTimeStep, maxCo, maxDeltaT, CoNum )
        
        runTime.increment()
        
        ref.ext_Info() << "Time = " << runTime.timeName() << ref.nl << ref.nl
        phi << aphiv_pos * rho_pos + aphiv_neg * rho_neg

        phiUp = ( aphiv_pos * rhoU_pos + aphiv_neg * rhoU_neg) + ( a_pos * p_pos + a_neg * p_neg ) * mesh.Sf()

        phiEp = aphiv_pos * ( rho_pos * ( e_pos + 0.5*U_pos.magSqr() ) + p_pos ) + aphiv_neg * ( rho_neg * ( e_neg + 0.5 * U_neg.magSqr() ) + p_neg )\
                + aSf * p_pos - aSf * p_neg
        
   
        muEff = turbulence.muEff()
        tauMC = ref.volTensorField( ref.word( "tauMC" ) , muEff * ref.fvc.grad(U).T().dev2() ) 

        # --- Solve density
        ref.solve( ref.fvm.ddt( rho ) + ref.fvc.div( phi ) )
        
        # --- Solve momentum
        ref.solve( ref.fvm.ddt( rhoU ) + ref.fvc.div( phiUp ) )
        
        U.dimensionedInternalField() << rhoU.dimensionedInternalField() / rho.dimensionedInternalField()
        U.correctBoundaryConditions()
        
        rhoU.ext_boundaryField() << rho.ext_boundaryField() * U.ext_boundaryField()
        
        rhoBydt = rho / runTime.deltaT()
        
        if not inviscid:
           solve( fvm.ddt( rho, U ) - fvc.ddt( rho, U ) - fvm.laplacian( muEff, U ) - fvc.div( tauMC ) )
           rhoU << rho * U
           pass
        
        # --- Solve energy
        sigmaDotU = ( ref.fvc.interpolate( muEff ) * mesh.magSf() * ref.fvc.snGrad( U ) + 
                      ( mesh.Sf() & ref.fvc.interpolate( tauMC ) ) ) & ( a_pos * U_pos + a_neg * U_neg )

        ref.solve( ref.fvm.ddt( rhoE ) + ref.fvc.div( phiEp ) - ref.fvc.div( sigmaDotU ) )
        
        e << rhoE() / rho() - 0.5 * U.magSqr() # mixed calculations
        e.correctBoundaryConditions()
        thermo.correct()

        rhoE.ext_boundaryField() << rho.ext_boundaryField() * ( e.ext_boundaryField() + 0.5 * U.ext_boundaryField().magSqr() )
        
        if not inviscid :
           k = man.volScalarField( ref.word( "k" ) , thermo.Cp() * muEff / Pr )

           # The initial C++ expression does not work properly, because of
           #  1. the order of expression arguments computation differs with C++
           #solve( fvm.ddt( rho, e ) - fvc.ddt( rho, e ) - fvm.laplacian( thermo.alpha(), e ) \
           #                                             + fvc.laplacian( thermo.alpha(), e ) - fvc.laplacian( k, T ) )

           solve( -fvc.laplacian( k, T ) + ( fvc.laplacian( turbulence.alpha(), e ) \
                                         + (- fvm.laplacian( turbulence.alphaEff(), e ) + (- fvc.ddt( rho, e ) + fvm.ddt( rho, e ) ) ) ) )
           
           thermo.correct()
           rhoE << rho * ( e + 0.5 * U.magSqr() )
           pass
        
        p.dimensionedInternalField() << rho.dimensionedInternalField() / psi.dimensionedInternalField()
        p.correctBoundaryConditions()

        rho.ext_boundaryField() << psi.ext_boundaryField() * p.ext_boundaryField() 
        
        turbulence.correct()
        runTime.write()

        ref.ext_Info() << "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << \
              "  ClockTime = " << runTime.elapsedClockTime() << " s" << ref.nl << ref.nl
        
        pass

    ref.ext_Info() << "End\n"

    import os
    return os.EX_OK


#--------------------------------------------------------------------------------------
import sys, os
from Foam import FOAM_VERSION
if FOAM_VERSION( ">=", "020100" ):
   if __name__ == "__main__" :
      argv = sys.argv
      os._exit( main_standalone( len( argv ), argv ) )
      pass
   pass   
else:
   from Foam.OpenFOAM import ext_Info
   ref.ext_Info()<< "\nTo use this solver, It is necessary to SWIG OpenFoam2.1.0 or higher \n "

   
#--------------------------------------------------------------------------------------

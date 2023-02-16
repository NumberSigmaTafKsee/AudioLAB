/*

    Cascaded IIR stages, header
    Copyright (C) 2002-2005 Jussi Laako

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/


#include <dsp/DSPOp.hh>


#ifndef IIRCASCADE_HH
    #define IIRCASCADE_HH


    /**
        Class to handle cascaded IIR stages as one filter.
    */
    class clIIRCascade
    {
            bool bInitialized;
        protected:
            long lStages;
            clDSPOp *IIR;
        public:
            clIIRCascade ();
            virtual ~clIIRCascade ();
            /**
                Initialize filter.
                
                \note You can re-initialize without uninitializing first.
                
                \param Coeffs Array of arrays of 5 coefficient bi-quad stages
                \param StageCount Number of cascaded bi-quad filter stages
            */
            bool Initialize (const float[][5], long);
            /// \overload
            bool Initialize (const double[][5], long);
            /**
                Uninitialize filter.
            */
            void Uninitialize ();
            /**
                Filter data in-place.
                
                \param Vect Source & destination vector
                \param Count Vector length
            */
            void Process (float *, long);
            /// \overload
            void Process (double *, long);
            /**
                Filter data out-of-place.
                
                \param Dest Destination vector
                \param Src Source vector
                \param Count Vector length
            */
            void Process (float *, const float *, long);
            void Process (double *, const double *, long);
            /**
                Clear filter feedback chain.
            */
            void Clear ();
    };


#endif

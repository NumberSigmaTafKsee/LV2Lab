#pragma once

// I dont know what happened here
// MDFM is in Amplifier
// This was some experiment I never finished.


#include "Bezier.hpp"
#include "DCFilter.hpp"
namespace Analog
{
    struct BezierDistortion
    {
        Bezier::Bezier<3> *cbezier;

        BezierDistortion() 
        {
            cbezier = nullptr;
            init();
        }
        ~BezierDistortion() {
            delete cbezier;
        }
        void init() {
            if(cbezier) delete cbezier;
            Random noise;
            int x1 = noise.randint(250,250);
            int y1 = noise.randint(0,250);
            int x2 = noise.randint(x1,500);
            int y2 = noise.randint(0,250);
            int x3 = noise.randint(x2,750);
            int y3 = noise.randint(0,250);
            int x4 = noise.randint(x3,1000);
            int y4 = noise.randint(0,250);
            cbezier = new Bezier::Bezier<3>(
            {                                
                {x1,y1},
                {x2,y2},
                {x3,y3},
                {x4,y4},
                
            });
        }
        DspFloatType Tick(DspFloatType I)
        {
            DspFloatType v = cbezier->valueAt(fabs(I),0)/1000.0;
            v = cbezier->valueAt(v,1)/250.0f;
            v *= signum(I);
            if(v < -1) v = -1.0f;
            if(v >  1) v = 1.0f;
            return v;
        }
    };
}
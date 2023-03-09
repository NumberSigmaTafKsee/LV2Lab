// generated automatically
// DO NOT MODIFY!
declare id "hornet";
declare name "Hornet";
declare category "Distortion";
declare description "Hornet simulation";

import("stdfaust.lib");

process =  fi.iir((b0/a0,b1/a0,b2/a0,b3/a0,b4/a0),(a1/a0,a2/a0,a3/a0,a4/a0)) : clip   with {
    LogPot(a, x) = ba.if(a, (exp(a * x) - 1) / (exp(a) - 1), x);
    Inverted(b, x) = ba.if(b, 1 - x, x);
    s = 0.993;
    fs = float(ma.SR);
    pre = _;
    a = 1.2715 - Fuzz ;
    clip(x) = (0.4 * (min(0.7514,max(-0.4514,x))));

        Volume = vslider("Volume[name:Volume]", 0.5, 0, 1, 0.01) : Inverted(0) : si.smooth(s);
    
        Sustain = vslider("Sustain[name:Sustain]", 0.5, 0, 1, 0.01) : Inverted(0) : si.smooth(s);
    
        Fuzz = vslider("Fuzz[name:Fuzz]", 0.5, 0, 1, 0.01) : Inverted(0) : si.smooth(s);
    
    b0 = Fuzz*(Fuzz*Volume*pow(fs,3)*(-3.36831187151837e-20*fs - 1.75582214579149e-16) + Volume*pow(fs,3)*(9.4649563589667e-21*fs + 4.93386022967413e-17)) + Sustain*(Fuzz*(Fuzz*Volume*pow(fs,3)*(2.89798007739403e-18*fs + 1.51064918927987e-14) + Volume*pow(fs,3)*(2.28655633153439e-18*fs + 1.19192830048069e-14)) + Volume*pow(fs,2)*(fs*(1.45891323583538e-19*fs + 7.40489463150802e-16) - 1.04296301457845e-13)) + Volume*pow(fs,2)*(fs*(2.4218162356217e-20*fs + 1.26859583357635e-16) + 3.21091305171869e-15);

    b1 = Fuzz*(Fuzz*Volume*pow(fs,3)*(1.34732474860735e-19*fs + 3.51164429158298e-16) + Volume*pow(fs,3)*(-3.78598254358668e-20*fs - 9.86772045934826e-17)) + Sustain*(Fuzz*(Fuzz*Volume*pow(fs,3)*(-1.15919203095761e-17*fs - 3.02129837855973e-14) + Volume*pow(fs,3)*(-9.14622532613756e-18*fs - 2.38385660096139e-14)) + Volume*pow(fs,3)*(-5.83565294334152e-19*fs - 1.4809789263016e-15)) + Volume*pow(fs,3)*(-9.6872649424868e-20*fs - 2.53719166715271e-16);

    b2 = Fuzz*(-2.02098712291102e-19*Fuzz*Volume*pow(fs,4) + 5.67897381538002e-20*Volume*pow(fs,4)) + Sustain*(Fuzz*(1.73878804643642e-17*Fuzz*Volume*pow(fs,4) + 1.37193379892063e-17*Volume*pow(fs,4)) + Volume*pow(fs,2)*(8.75347941501228e-19*pow(fs,2) + 2.08592602915691e-13)) + Volume*pow(fs,2)*(1.45308974137302e-19*pow(fs,2) - 6.42182610343738e-15);

    b3 = Fuzz*(Fuzz*Volume*pow(fs,3)*(1.34732474860735e-19*fs - 3.51164429158298e-16) + Volume*pow(fs,3)*(-3.78598254358668e-20*fs + 9.86772045934826e-17)) + Sustain*(Fuzz*(Fuzz*Volume*pow(fs,3)*(-1.15919203095761e-17*fs + 3.02129837855973e-14) + Volume*pow(fs,3)*(-9.14622532613756e-18*fs + 2.38385660096139e-14)) + Volume*pow(fs,3)*(-5.83565294334152e-19*fs + 1.4809789263016e-15)) + Volume*pow(fs,3)*(-9.6872649424868e-20*fs + 2.53719166715271e-16);

    b4 = Fuzz*(Fuzz*Volume*pow(fs,3)*(-3.36831187151837e-20*fs + 1.75582214579149e-16) + Volume*pow(fs,3)*(9.4649563589667e-21*fs - 4.93386022967413e-17)) + Sustain*(Fuzz*(Fuzz*Volume*pow(fs,3)*(2.89798007739403e-18*fs - 1.51064918927987e-14) + Volume*pow(fs,3)*(2.28655633153439e-18*fs - 1.19192830048069e-14)) + Volume*pow(fs,2)*(fs*(1.45891323583538e-19*fs - 7.40489463150802e-16) - 1.04296301457845e-13)) + Volume*pow(fs,2)*(fs*(2.4218162356217e-20*fs - 1.26859583357635e-16) + 3.21091305171869e-15);

    a0 = Fuzz*(Fuzz*fs*(fs*(fs*(-2.57087433571955e-21*fs - 3.20282580029198e-16) - 1.59955479510613e-12) - 1.59265781983301e-11) + fs*(fs*(fs*(7.22415688337201e-22*fs + 1.15708148345401e-16) + 5.83730418035165e-13) + 5.81517816570128e-12)) + Sustain*(Fuzz*(Fuzz*fs*(fs*(fs*(-2.22731835703847e-20*fs - 1.15680872556631e-16) - 1.02110226030461e-14) - 8.86189120121937e-14) + fs*(fs*(fs*(2.63346198155234e-19*fs + 1.39560593212619e-15) + 1.23553821745518e-13) + 1.09557219057811e-12)) + fs*(fs*(fs*(2.00860283725342e-19*fs + 1.06638545237487e-15) + 1.10730441129884e-13) + 2.30413820563986e-12) + 1.30107041069324e-11) + fs*(fs*(fs*(1.84845864738235e-21*fs + 2.48814775695488e-16) + 1.25274019677064e-12) + 4.3094216014379e-11) + 3.05506357605318e-10;

    a1 = Fuzz*(Fuzz*fs*(pow(fs,2)*(1.02834973428782e-20*fs + 6.40565160058397e-16) - 3.18531563966602e-11) + fs*(pow(fs,2)*(-2.8896627533488e-21*fs - 2.31416296690802e-16) + 1.16303563314026e-11)) + Sustain*(Fuzz*(Fuzz*fs*(pow(fs,2)*(8.9092734281539e-20*fs + 2.31361745113263e-16) - 1.77237824024387e-13) + fs*(pow(fs,2)*(-1.05338479262093e-18*fs - 2.79121186425238e-15) + 2.19114438115622e-12)) + fs*(pow(fs,2)*(-8.03441134901368e-19*fs - 2.13277090474974e-15) + 4.60827641127972e-12) + 5.20428164277295e-11) + fs*(pow(fs,2)*(-7.39383458952941e-21*fs - 4.97629551390976e-16) + 8.6188432028758e-11) + 1.22202543042127e-9;

    a2 = Fuzz*(Fuzz*pow(fs,2)*(-1.54252460143173e-20*pow(fs,2) + 3.19910959021226e-12) + pow(fs,2)*(4.33449413002321e-21*pow(fs,2) - 1.16746083607033e-12)) + Sustain*(Fuzz*(Fuzz*pow(fs,2)*(-1.33639101422308e-19*pow(fs,2) + 2.04220452060922e-14) + pow(fs,2)*(1.5800771889314e-18*pow(fs,2) - 2.47107643491035e-13)) + pow(fs,2)*(1.20516170235205e-18*pow(fs,2) - 2.21460882259768e-13) + 7.80642246415943e-11) + pow(fs,2)*(1.10907518842941e-20*pow(fs,2) - 2.50548039354128e-12) + 1.83303814563191e-9;

    a3 = Fuzz*(Fuzz*fs*(pow(fs,2)*(1.02834973428782e-20*fs - 6.40565160058397e-16) + 3.18531563966602e-11) + fs*(pow(fs,2)*(-2.8896627533488e-21*fs + 2.31416296690802e-16) - 1.16303563314026e-11)) + Sustain*(Fuzz*(Fuzz*fs*(pow(fs,2)*(8.9092734281539e-20*fs - 2.31361745113263e-16) + 1.77237824024387e-13) + fs*(pow(fs,2)*(-1.05338479262093e-18*fs + 2.79121186425238e-15) - 2.19114438115622e-12)) + fs*(pow(fs,2)*(-8.03441134901368e-19*fs + 2.13277090474974e-15) - 4.60827641127972e-12) + 5.20428164277295e-11) + fs*(pow(fs,2)*(-7.39383458952941e-21*fs + 4.97629551390976e-16) - 8.6188432028758e-11) + 1.22202543042127e-9;

    a4 = Fuzz*(Fuzz*fs*(fs*(fs*(-2.57087433571955e-21*fs + 3.20282580029198e-16) - 1.59955479510613e-12) + 1.59265781983301e-11) + fs*(fs*(fs*(7.22415688337201e-22*fs - 1.15708148345401e-16) + 5.83730418035165e-13) - 5.81517816570128e-12)) + Sustain*(Fuzz*(Fuzz*fs*(fs*(fs*(-2.22731835703847e-20*fs + 1.15680872556631e-16) - 1.02110226030461e-14) + 8.86189120121937e-14) + fs*(fs*(fs*(2.63346198155234e-19*fs - 1.39560593212619e-15) + 1.23553821745518e-13) - 1.09557219057811e-12)) + fs*(fs*(fs*(2.00860283725342e-19*fs - 1.06638545237487e-15) + 1.10730441129884e-13) - 2.30413820563986e-12) + 1.30107041069324e-11) + fs*(fs*(fs*(1.84845864738235e-21*fs - 2.48814775695488e-16) + 1.25274019677064e-12) - 4.3094216014379e-11) + 3.05506357605318e-10;
};
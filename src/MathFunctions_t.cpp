/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "MathFunctions.hpp"
#include "RandomGenerator.hpp"
#include "catch2/catch.hpp"
#include <algorithm>

TEST_CASE("CMathFunctions::calcWVMie::BH", "[MathFunctions][CMathFunctions]")
{
    // compare results with
    // Bohren & Huffman (1998), Absorption and Scattering of Light by Small Particles
    // see Appendix A, page 482

    uint nr_theta = 21;
    dcomplex ri {1.55, 0.0};
    double x {PIx2 * 0.525 / 0.6328};
    double qext, qabs, qsca, gsca {};
    double S11[nr_theta], S12[nr_theta], S33[nr_theta], S34[nr_theta];

    dlist theta(nr_theta);
    for(uint th = 0; th < nr_theta; th++) {
        theta[th] = th * PI / double(nr_theta - 1);
    }

    bool res {CMathFunctions::calcWVMie(x, theta, ri, qext, qabs, qsca, gsca, S11, S12, S33, S34)};

    REQUIRE(res);

    REQUIRE(std::abs(qext - 3.10543) < 1e-5);
    REQUIRE(std::abs(qabs - 0.00000) < 1e-5);
    REQUIRE(std::abs(qsca - 3.10543) < 1e-5);
    // REQUIRE(abs(qbk - 2.92534) < 1e-5);

    double BH_S11[] = {
        0.100000e+01, 0.785390e+00, 0.356897e+00, 0.766119e-01, 0.355355e-01, 0.701845e-01, 0.574313e-01, 0.219660e-01, 0.125959e-01, 0.173750e-01,
        0.124601e-01, 0.679093e-02, 0.954239e-02, 0.863419e-02, 0.227421e-02, 0.543998e-02, 0.160243e-01, 0.188852e-01, 0.195254e-01, 0.301676e-01,
        0.383189e-01};
    double BH_S12[] = {
        0.000000e+00,-0.459811e-02,-0.458541e-01,-0.364744e+00,-0.534997e+00, 0.959953e-02, 0.477927e-01,-0.440604e+00,-0.831996e+00, 0.341670e-01,
        0.230462e+00,-0.713472e+00,-0.756255e+00,-0.281215e+00,-0.239612e+00,-0.850804e+00,-0.706334e+00,-0.891081e+00,-0.783319e+00,-0.196194e+00,
        0.000000e+00};
    double BH_S33[] = {
        0.100000e+01, 0.999400e+00, 0.986022e+00, 0.843603e+00, 0.686967e+00, 0.959825e+00, 0.985371e+00, 0.648043e+00, 0.203255e+00, 0.795354e+00,
        0.937497e+00,-0.717397e-02,-0.394748e-01, 0.536251e+00, 0.967602e+00, 0.187531e+00, 0.495254e+00, 0.453277e+00,-0.391613e+00,-0.962069e+00,
       -0.100000e+01};
    double BH_S34[] = {
        0.000000e+00, 0.343261e-01, 0.160184e+00, 0.394076e+00,-0.491787e+00,-0.280434e+00, 0.163584e+00, 0.621216e+00,-0.516208e+00,-0.605182e+00,
        0.260742e+00, 0.700647e+00,-0.653085e+00,-0.795835e+00, 0.795798e-01,-0.490882e+00,-0.505781e+00,-0.226817e-01, 0.482752e+00, 0.189556e+00,
        0.000000e+00};

    for(uint th = 0; th < nr_theta; th++) {
        REQUIRE(std::abs( S11[th] / S11[0]  - BH_S11[th]) < 1e-6);
        REQUIRE(std::abs(-S12[th] / S11[th] - BH_S12[th]) < 1e-6);
        REQUIRE(std::abs( S33[th] / S11[th] - BH_S33[th]) < 1e-6);
        REQUIRE(std::abs( S34[th] / S11[th] - BH_S34[th]) < 1e-6);
    }
}

TEST_CASE("CMathFunctions::calcWVMie::W1", "[MathFunctions][CMathFunctions]")
{
    // compare results with
    // Wiscombe (1979), Mie Scattering Calculations: Advances in Technique and Fast, Vector-speed Computer Codes
    // see Appendix, A22

    uint nr_theta = 37;
    dcomplex ri {1.5, 0.0};
    auto x = GENERATE(10.0, 100.0, 1000.0, 5000.0);
    double qext, qabs, qsca, gsca {};
    double S11[nr_theta], S12[nr_theta], S33[nr_theta], S34[nr_theta];

    dlist theta(nr_theta);
    for(uint th = 0; th < nr_theta; th++) {
        theta[th] = th * PI / double(nr_theta - 1);
    }

    bool res {CMathFunctions::calcWVMie(x, theta, ri, qext, qabs, qsca, gsca, S11, S12, S33, S34)};
    CAPTURE(x);
    REQUIRE(res);

    if(x == 10.0) {
        REQUIRE(std::abs(qext - 2.881999) < 1e-6);
        REQUIRE(std::abs(qabs - 0.000000) < 1e-6);
        REQUIRE(std::abs(qsca - 2.881999) < 1e-6);
        REQUIRE(std::abs(gsca - 0.742913) < 1e-6);

        double W_S11[] = {
            0.520856e+04, 0.428132e+04, 0.234775e+04, 0.818121e+03, 0.149547e+03, 0.139195e+02, 0.768071e+02, 0.128393e+03, 0.789487e+02, 0.293697e+02,
            0.518546e+02, 0.631908e+02, 0.341567e+02, 0.256714e+02, 0.320167e+02, 0.235919e+02, 0.194168e+02, 0.184988e+02, 0.917521e+01, 0.963035e+01,
            0.136986e+02, 0.485859e+01, 0.280440e+01, 0.853248e+01, 0.439826e+01, 0.570995e+01, 0.113872e+02, 0.304142e+01, 0.675713e+01, 0.242055e+02,
            0.159589e+02, 0.161562e+02, 0.582167e+02, 0.679231e+02, 0.286866e+02, 0.249122e+02, 0.423766e+02};
        double W_S12[] = {
            0.0000,-0.0239,-0.0963,-0.2080,-0.2887, 0.4685, 0.0005,-0.0892,-0.0767, 0.7251,
            0.4834,-0.0592,-0.0163, 0.9793, 0.3242,-0.3951, 0.5014, 0.9514,-0.0269, 0.1127,
            0.7525, 0.7257, 0.8993, 0.5987,-0.4844, 0.3364, 0.9396, 0.6460, 0.6282, 0.9286,
            0.7664,-0.4108, 0.2385, 0.5276, 0.6986, 0.0830, 0.0000};

        for(uint th = 0; th < nr_theta; th++) {
            REQUIRE(std::abs(1.0 - S11[th] / W_S11[th]) < 1e-5);
            REQUIRE(std::abs(S12[th] / S11[th] - W_S12[th]) < 1e-4);
        }
    } else if(x == 100.0) {
        REQUIRE(std::abs(qext - 2.094388) < 1e-6);
        REQUIRE(std::abs(qabs - 0.000000) < 1e-6);
        REQUIRE(std::abs(qsca - 2.094388) < 1e-6);
        REQUIRE(std::abs(gsca - 0.818246) < 1e-6);

        double W_S11[] = {
            0.275508e+08, 0.106440e+06, 0.586911e+05, 0.258367e+05, 0.202131e+05, 0.192755e+05, 0.159327e+05, 0.538269e+04, 0.866710e+04, 0.540267e+04,
            0.305275e+04, 0.510959e+04, 0.167021e+04, 0.202680e+04, 0.960828e+03, 0.122361e+04, 0.394465e+03, 0.254063e+03, 0.505356e+03, 0.507283e+03,
            0.131614e+03, 0.222694e+03, 0.133922e+03, 0.705029e+02, 0.401041e+02, 0.100503e+03, 0.167811e+03, 0.232785e+03, 0.243533e+03, 0.235398e+03,
            0.186086e+03, 0.121783e+04, 0.353283e+04, 0.239888e+04, 0.272592e+04, 0.622951e+03, 0.434048e+04};
        double W_S12[] = {
            0.0000,-0.0107, 0.0098, 0.0712,-0.0248,-0.0730,-0.0540, 0.2001, 0.0919,-0.1301,
            0.4545, 0.1026, 0.4679,-0.2130, 0.7357,-0.0054, 0.8365,-0.2535,-0.0016,-0.4214,
           -0.5649,-0.3412, 0.3830,-0.0700,-0.8726, 0.1822, 0.5235, 0.4015, 0.3292, 0.0578,
            0.5225,-0.8534,-0.8664, 0.6181, 0.8582, 0.5069, 0.0000};

        for(uint th = 0; th < nr_theta; th++) {
            REQUIRE(std::abs(1.0 - S11[th] / W_S11[th]) < 1e-5);
            REQUIRE(std::abs(S12[th] / S11[th] - W_S12[th]) < 1e-4);
        }
    } else if(x == 1000.0) {
        REQUIRE(std::abs(qext - 2.013945) < 1e-6);
        REQUIRE(std::abs(qabs - 0.000000) < 1e-6);
        REQUIRE(std::abs(qsca - 2.013945) < 1e-6);
        REQUIRE(std::abs(gsca - 0.827882) < 1e-6);

        double W_S11[] = {
            0.253568e+12, 0.482785e+07, 0.883654e+06, 0.163897e+07, 0.136169e+07, 0.113470e+07, 0.103882e+07, 0.929666e+06, 0.863859e+06, 0.766220e+06,
            0.266763e+06, 0.324223e+06, 0.282935e+06, 0.237449e+06, 0.160593e+06, 0.525112e+05, 0.789951e+05, 0.306843e+05, 0.248427e+05, 0.214530e+05,
            0.106283e+05, 0.152392e+05, 0.197453e+05, 0.912219e+04, 0.714740e+04, 0.791914e+04, 0.952024e+04, 0.109979e+05, 0.666722e+04, 0.168603e+05,
            0.211737e+05, 0.559288e+04, 0.208814e+06, 0.267605e+06, 0.392870e+06, 0.224420e+06, 0.257577e+07};
        double W_S12[] = {
            0.0000,-0.0455,-0.1300,-0.0275, 0.0909, 0.0469, 0.0701, 0.0406,-0.0202,-0.0694,
            0.5051, 0.1317,-0.0010,-0.1018,-0.0523, 0.3454,-0.0192,-0.3916, 0.0127,-0.8087,
           -0.8588,-0.6411,-0.6625,-0.5131, 0.3474,-0.4821,-0.2850,-0.2910, 0.3604,-0.3089,
           -0.4512, 0.4230,-0.0666, 0.0290,-0.8299,-0.0267, 0.0000};

        for(uint th = 0; th < nr_theta; th++) {
            REQUIRE(std::abs(1.0 - S11[th] / W_S11[th]) < 1e-5);
            REQUIRE(std::abs(S12[th] / S11[th] - W_S12[th]) < 1e-4);
        }
    } else if(x == 5000.0) {
        REQUIRE(std::abs(qext - 2.008650) < 1e-6);
        REQUIRE(std::abs(qabs - 0.000000) < 1e-6);
        REQUIRE(std::abs(qsca - 2.008650) < 1e-6);
        REQUIRE(std::abs(gsca - 0.829592) < 1e-6);

        double W_S11[] = {
            0.157609e+15, 0.394653e+08, 0.229931e+08, 0.339037e+08, 0.262577e+08, 0.236482e+08, 0.155814e+08, 0.323654e+08, 0.234225e+08, 0.178268e+08,
            0.673414e+07, 0.586469e+07, 0.709207e+07, 0.262123e+07, 0.310573e+07, 0.158540e+07, 0.219096e+07, 0.101618e+07, 0.567675e+06, 0.323205e+06,
            0.241279e+06, 0.249034e+06, 0.143024e+06, 0.417106e+06, 0.253470e+06, 0.199198e+06, 0.273150e+06, 0.254040e+06, 0.185353e+06, 0.363299e+06,
            0.272294e+06, 0.256060e+06, 0.130619e+08, 0.372204e+07, 0.453092e+07, 0.205242e+07, 0.237786e+09};
        double W_S12[] = {
            0.0000,-0.0223, 0.1179, 0.1247, 0.0279, 0.1168, 0.2597,-0.0712,-0.0175,-0.0512,
            0.5594, 0.3873, 0.0148, 0.8293, 0.0368, 0.0580,-0.1909,-0.3772,-0.4089,-0.7399,
           -0.7719,-0.6412, 0.9228,-0.8142,-0.1234, 0.0435,-0.3815,-0.1249,-0.1602,-0.6980,
           -0.7087,-0.0846,-0.8202,-0.1764, 0.2883, 0.6997, 0.0000};

        for(uint th = 0; th < nr_theta; th++) {
            REQUIRE(std::abs(1.0 - S11[th] / W_S11[th]) < 1e-5);
            REQUIRE(std::abs(S12[th] / S11[th] - W_S12[th]) < 1e-4);
        }
    }
}

TEST_CASE("CMathFunctions::calcWVMie::W2", "[MathFunctions][CMathFunctions]")
{
    // compare results with
    // Wiscombe (1979), Mie Scattering Calculations: Advances in Technique and Fast, Vector-speed Computer Codes
    // see Appendix, A22

    uint nr_theta = 37;
    dcomplex ri {1.5, 0.1};
    auto x = GENERATE(10.0, 100.0, 1000.0, 5000.0);
    double qext, qabs, qsca, gsca {};
    double S11[nr_theta], S12[nr_theta], S33[nr_theta], S34[nr_theta];

    dlist theta(nr_theta);
    for(uint th = 0; th < nr_theta; th++) {
        theta[th] = th * PI / double(nr_theta - 1);
    }

    bool res {CMathFunctions::calcWVMie(x, theta, ri, qext, qabs, qsca, gsca, S11, S12, S33, S34)};
    CAPTURE(x);
    REQUIRE(res);

    if(x == 10.0) {
        REQUIRE(std::abs(qext - 2.459791) < 1e-6);
        REQUIRE(std::abs(qabs - 1.224646) < 1e-6);
        REQUIRE(std::abs(qsca - 1.235144) < 1e-6);
        REQUIRE(std::abs(gsca - 0.922350) < 1e-6);

        double W_S11[] = {
            0.379171e+04, 0.300320e+04, 0.141624e+04, 0.313014e+03, 0.124235e+02, 0.296988e+02, 0.273164e+02, 0.112113e+02, 0.109517e+02, 0.607843e+01,
            0.220902e+01, 0.632075e+01, 0.646946e+01, 0.225394e+01, 0.215826e+01, 0.392848e+01, 0.299433e+01, 0.163623e+01, 0.183556e+01, 0.209544e+01,
            0.166228e+01, 0.137914e+01, 0.153058e+01, 0.147431e+01, 0.116521e+01, 0.135300e+01, 0.174359e+01, 0.136826e+01, 0.798073e+00, 0.974236e+00,
            0.133396e+01, 0.141816e+01, 0.148012e+01, 0.126487e+01, 0.106733e+01, 0.172292e+01, 0.231818e+01};
        double W_S12[] = {
            0.0000,-0.0014,-0.0068,-0.0301,-0.6183,-0.2518,-0.2817,-0.7359,-0.6378,-0.5132,
           -0.9235,-0.7194,-0.6077,-0.1744,-0.2426,-0.7454,-0.6373, 0.3019,-0.0893,-0.8614,
           -0.6653, 0.2706, 0.0790,-0.7132,-0.8966, 0.1033, 0.3819,-0.0370,-0.6271, 0.0599,
            0.3753, 0.1218,-0.2643,-0.6463,-0.8175,-0.2177, 0.0000};

        for(uint th = 0; th < nr_theta; th++) {
            REQUIRE(std::abs(1.0 - S11[th] / W_S11[th]) < 1e-5);
            REQUIRE(std::abs(S12[th] / S11[th] - W_S12[th]) < 1e-4);
        }
    } else if(x == 100.0) {
        REQUIRE(std::abs(qext - 2.089822) < 1e-6);
        REQUIRE(std::abs(qabs - 0.957688) < 1e-6);
        REQUIRE(std::abs(qsca - 1.132134) < 1e-6);
        REQUIRE(std::abs(gsca - 0.950392) < 1e-6);

        double W_S11[] = {
            0.273645e+08, 0.911574e+05, 0.130172e+05, 0.314304e+04, 0.121824e+04, 0.911319e+03, 0.801673e+03, 0.629347e+03, 0.465786e+03, 0.370932e+03,
            0.317391e+03, 0.269848e+03, 0.230451e+03, 0.202758e+03, 0.180401e+03, 0.162813e+03, 0.149203e+03, 0.138475e+03, 0.130113e+03, 0.123599e+03,
            0.118559e+03, 0.114659e+03, 0.111687e+03, 0.109437e+03, 0.107749e+03, 0.106503e+03, 0.105600e+03, 0.104961e+03, 0.104521e+03, 0.104229e+03, 
            0.104044e+03, 0.103935e+03, 0.103877e+03, 0.103850e+03, 0.103840e+03, 0.103837e+03, 0.103837e+03};
        double W_S12[] = {
            0.0000,-0.0247,-0.0848,-0.2375,-0.4798,-0.5484,-0.5537,-0.6268,-0.7490,-0.8418,
           -0.8905,-0.9395,-0.9797,-0.9960,-0.9944,-0.9750,-0.9395,-0.8885,-0.8273,-0.7576,
           -0.6831,-0.6072,-0.5316,-0.4586,-0.3897,-0.3257,-0.2673,-0.2148,-0.1684,-0.1278,
           -0.0932,-0.0643,-0.0409,-0.0229,-0.0101,-0.0025, 0.0000};

        for(uint th = 0; th < nr_theta; th++) {
            REQUIRE(std::abs(1.0 - S11[th] / W_S11[th]) < 1e-5);
            REQUIRE(std::abs(S12[th] / S11[th] - W_S12[th]) < 1e-4);
        }
    } else if(x == 1000.0) {
        REQUIRE(std::abs(qext - 2.019703) < 1e-6);
        REQUIRE(std::abs(qabs - 0.912770) < 1e-6);
        REQUIRE(std::abs(qsca - 1.106932) < 1e-6);
        REQUIRE(std::abs(gsca - 0.950880) < 1e-6);

        double W_S11[] = {
            0.255002e+12, 0.103489e+07, 0.270449e+06, 0.145265e+06, 0.101955e+06, 0.803929e+05, 0.646676e+05, 0.528211e+05, 0.436516e+05, 0.364909e+05,
            0.308618e+05, 0.264252e+05, 0.229245e+05, 0.201609e+05, 0.179799e+05, 0.162601e+05, 0.149061e+05, 0.138426e+05, 0.130097e+05, 0.123601e+05,
            0.118559e+05, 0.114671e+05, 0.111696e+05, 0.109441e+05, 0.107752e+05, 0.106505e+05, 0.105601e+05, 0.104961e+05, 0.104520e+05, 0.104227e+05,
            0.104042e+05, 0.103933e+05, 0.103874e+05, 0.103846e+05, 0.103836e+05, 0.103834e+05, 0.103834e+05};
        double W_S12[] = {
            0.0000,-0.0682,-0.1681,-0.2830,-0.3931,-0.4861,-0.5789,-0.6682,-0.7517,-0.8267,
           -0.8909,-0.9417,-0.9765,-0.9939,-0.9928,-0.9737,-0.9379,-0.8878,-0.8265,-0.7571,
           -0.6829,-0.6069,-0.5315,-0.4586,-0.3897,-0.3258,-0.2674,-0.2149,-0.1684,-0.1279,
           -0.0932,-0.0643,-0.0409,-0.0229,-0.0101,-0.0025, 0.0000};

        for(uint th = 0; th < nr_theta; th++) {
            REQUIRE(std::abs(1.0 - S11[th] / W_S11[th]) < 1e-5);
            REQUIRE(std::abs(S12[th] / S11[th] - W_S12[th]) < 1e-4);
        }
    } else if(x == 5000.0) {
        REQUIRE(std::abs(qext - 2.006775) < 1e-6);
        REQUIRE(std::abs(qabs - 0.907582) < 1e-6);
        REQUIRE(std::abs(qsca - 1.099193) < 1e-6);
        REQUIRE(std::abs(gsca - 0.950650) < 1e-6);

        double W_S11[] = {
            0.157315e+15, 0.772728e+07, 0.417917e+07, 0.311291e+07, 0.245545e+07, 0.197572e+07, 0.160555e+07, 0.131668e+07, 0.109011e+07, 0.911909e+06,
            0.771474e+06, 0.660667e+06, 0.573177e+06, 0.504089e+06, 0.449555e+06, 0.406550e+06, 0.372691e+06, 0.346095e+06, 0.325266e+06, 0.309020e+06,
            0.296412e+06, 0.286689e+06, 0.279248e+06, 0.273608e+06, 0.269384e+06, 0.266266e+06, 0.264005e+06, 0.262403e+06, 0.261300e+06, 0.260568e+06,
            0.260105e+06, 0.259832e+06, 0.259684e+06, 0.259616e+06, 0.259591e+06, 0.259585e+06, 0.259585e+06};
        double W_S12[] = {
            0.0000,-0.1103,-0.1975,-0.2927,-0.3902,-0.4856,-0.5788,-0.6680,-0.7513,-0.8264,
           -0.8906,-0.9414,-0.9763,-0.9936,-0.9926,-0.9735,-0.9378,-0.8878,-0.8264,-0.7570,
           -0.6829,-0.6069,-0.5315,-0.4586,-0.3897,-0.3258,-0.2674,-0.2149,-0.1684,-0.1279,
           -0.0932,-0.0643,-0.0409,-0.0229,-0.0101,-0.0025, 0.0000};

        for(uint th = 0; th < nr_theta; th++) {
            REQUIRE(std::abs(1.0 - S11[th] / W_S11[th]) < 1e-5);
            REQUIRE(std::abs(S12[th] / S11[th] - W_S12[th]) < 1e-4);
        }
    }
}

TEST_CASE("CMathFunctions::findRootBrent::Phi", "[MathFunctions][CMathFunctions]")
{
    CRandomGenerator randgen {};
    randgen.init(123);

    auto phipar = GENERATE(-1.0, -0.5, 0.0, 0.5, 1.0);

    double phi {};
    double phi_mean {0.0};
    double phi_min {PIx2};
    double phi_max {0.0};
    constexpr std::size_t n_samples {1024 * 1024};

    for (uint i = 0; i < 512; i++) {
        double phi_rnd {randgen.getRND() * PIx2};
        double rnd {CMathFunctions::getPhiIntegral(phi_rnd, {phipar, 0.0})};
        phi = CMathFunctions::findRootBrent(0.0, PIx2, &CMathFunctions::getPhiIntegral, {phipar, rnd});
        // the procedure should find the root of phi
        REQUIRE(std::abs(phi_rnd - phi) < 1e-4);
    }

    for (std::size_t i = 0; i < n_samples; ++i) {
        phi = CMathFunctions::findRootBrent(0.0, PIx2, &CMathFunctions::getPhiIntegral, {phipar, randgen.getRND()});
        phi_mean += phi;
        phi_min = std::min(phi, phi_min);
        phi_max = std::max(phi, phi_max);
    }

    CAPTURE(phipar);

    phi_mean /= static_cast<double>(n_samples);
    // phi should be between 0 and 2pi
    REQUIRE(phi_min >= 0.0);
    REQUIRE(phi_max <= PIx2);
    REQUIRE(phi_min < phi_max);
    // mean value should be approximately pi
    REQUIRE(std::abs(phi_mean - PI) < 0.01);
}

TEST_CASE("CMathFunctions::findRootBrent::DHG", "[MathFunctions][CMathFunctions]")
{
    CRandomGenerator randgen {};
    randgen.init(123);

    auto g = GENERATE(-0.99, -0.5, 0.0, 0.5, 0.99);
    auto alpha = GENERATE(0.0, 0.33, 0.67, 1.0);

    // mu = cos(theta)
    double mu {};
    double mu_mean {0.0};
    double mu_min {1.0};
    double mu_max {-1.0};
    constexpr std::size_t n_samples {1024 * 1024};

    for (uint i = 0; i < 512; i++) {
        double mu_rnd {2.0 * randgen.getRND() - 1.0};
        double rnd {CMathFunctions::getDHGIntegral(mu_rnd, {g, alpha, 0.0})};
        mu = CMathFunctions::findRootBrent(-1.0, 1.0, &CMathFunctions::getDHGIntegral, {g, alpha, rnd});
        // the procedure should find the root of mu
        REQUIRE(std::abs(mu_rnd - mu) < 1e-4);
    }

    for (std::size_t i = 0; i < n_samples; ++i) {
        mu = CMathFunctions::findRootBrent(-1.0, 1.0, &CMathFunctions::getDHGIntegral, {g, alpha, randgen.getRND()});
        mu_mean += mu;
        mu_min = std::min(mu, mu_min);
        mu_max = std::max(mu, mu_max);
    }

    CAPTURE(g);
    CAPTURE(alpha);

    mu_mean /= static_cast<double>(n_samples);
    // mu should be between -1 and 1
    REQUIRE(mu_min >= -1.0);
    REQUIRE(mu_max <= 1.0);
    REQUIRE(mu_min < mu_max);
    // mean value should be approximately Eq. (A2) in Draine 2003, ApJ 598, 1017
    double first_moment {g * (1.0 + alpha * (3.0 + 2.0 * g * g) / 5.0) / (1.0 + alpha * (1.0 + 2.0 * g * g) / 3.0)};
    REQUIRE(std::abs(mu_mean - first_moment) < 0.01);
}

TEST_CASE("CMathFunctions::findRootBrent::TTHG", "[MathFunctions][CMathFunctions]")
{
    CRandomGenerator randgen {};
    randgen.init(123);

    auto g1 = GENERATE(0.0, 0.5, 0.99);
    auto g2 = GENERATE(-0.99, -0.5, 0.0);
    auto weight = GENERATE(0.0, 0.33, 0.67, 1.0);

    // mu = cos(theta)
    double mu {};
    double mu_mean {0.0};
    double mu_min {1.0};
    double mu_max {-1.0};
    constexpr std::size_t n_samples {1024 * 1024};

    for (uint i = 0; i < 1024; i++) {
        double mu_rnd {2.0 * randgen.getRND() - 1.0};
        double rnd {CMathFunctions::getTTHGIntegral(mu_rnd, {g1, g2, weight, 0.0})};
        mu = CMathFunctions::findRootBrent(-1.0, 1.0, &CMathFunctions::getTTHGIntegral, {g1, g2, weight, rnd});
        // the procedure should find the root of mu
        REQUIRE(std::abs(mu_rnd - mu) < 1e-4);
    }

    for (std::size_t i = 0; i < n_samples; ++i) {
        mu = CMathFunctions::findRootBrent(-1.0, 1.0, &CMathFunctions::getTTHGIntegral, {g1, g2, weight, randgen.getRND()});
        mu_mean += mu;
        mu_min = std::min(mu, mu_min);
        mu_max = std::max(mu, mu_max);
    }

    CAPTURE(g1);
    CAPTURE(g2);
    CAPTURE(weight);

    mu_mean /= static_cast<double>(n_samples);
    // mu should be between -1 and 1
    REQUIRE(mu_min >= -1.0);
    REQUIRE(mu_max <= 1.0);
    REQUIRE(mu_min < mu_max);
    // mean value should be approximately w*g1 + (1-w)*g2
    double first_moment {weight * g1 + (1.0 - weight) * g2};
    REQUIRE(std::abs(mu_mean - first_moment) < 0.01);
}

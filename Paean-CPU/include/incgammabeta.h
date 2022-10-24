#ifndef PAEAN_INCGAMMABETA_H
#define PAEAN_INCGAMMABETA_H

#include <math.h>
#include <stdio.h>
#include <cfloat>
#include <limits>

inline double SQR(double x) {
    return x * x;
}
int ngau = 18;
double gauy[18] = {
    0.0021695375159141994, 0.011413521097787704, 0.027972308950302116,
    0.051727015600492421, 0.082502225484340941, 0.12007019910960293,
    0.16415283300752470, 0.21442376986779355, 0.27051082840644336,
    0.33199876341447887, 0.39843234186401943, 0.46931971407375483,
    0.54413605556657973, 0.62232745288031077, 0.70331500465597174,
    0.78649910768313447, 0.87126389619061517, 0.95698180152629142};
double gauw[18] = {
    0.0055657196642445571, 0.012915947284065419, 0.020181515297735382,
    0.027298621498568734, 0.034213810770299537, 0.040875750923643261,
    0.047235083490265582, 0.053244713977759692, 0.058860144245324798,
    0.064039797355015485, 0.068745323835736408, 0.072941885005653087,
    0.076598410645870640, 0.079687828912071670, 0.082187266704339706,
    0.084078218979661945, 0.085346685739338721, 0.085983275670394821};

int SWITCH = 3000;
double EPS = DBL_EPSILON;
double FPMIN = DBL_MIN / DBL_EPSILON;

double betai(const double a, const double b, const double x);
static double betacf(const double a, const double b, const double x);
double betaiapprox(double a, double b, double x);
double invbetai(double p, double a, double b);

double betai(const double a, const double b, const double x) {
    double bt;
    //if (a <= 0.0 || b <= 0.0) throw("Bad a or b in routine betai");
    //if (x < 0.0 || x > 1.0) throw("Bad x in routine betai");
    if (x == 0.0 || x == 1.0)
        return x;
    if (a > SWITCH && b > SWITCH)
        return betaiapprox(a, b, x);
    bt = exp(lgammaf(a + b) - lgammaf(a) - lgammaf(b) + a * log(x) +
             b * log(1.0 - x));
    if (x < (a + 1.0) / (a + b + 2.0))
        return bt * betacf(a, b, x) / a;
    else
        return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
}

double betacf(const double a, const double b, const double x) {
    int m, m2;
    double aa, c, d, del, h, qab, qam, qap;
    qab = a + b;
    qap = a + 1.0;
    qam = a - 1.0;
    c = 1.0;
    d = 1.0 - qab * x / qap;
    if (fabs(d) < FPMIN)
        d = FPMIN;
    d = 1.0 / d;
    h = d;
    //#pragma unroll
    for (m = 1; m < 100; m++) {
        m2 = 2 * m;
        aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < FPMIN)
            d = FPMIN;
        c = 1.0 + aa / c;
        if (fabs(c) < FPMIN)
            c = FPMIN;
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < FPMIN)
            d = FPMIN;
        c = 1.0 + aa / c;
        if (fabs(c) < FPMIN)
            c = FPMIN;
        d = 1.0 / d;
        del = d * c;
        h *= del;
        if (fabs(del - 1.0) <= EPS)
            break;
    }
    return h;
}

double betaiapprox(double a, double b, double x) {
    int j;
    double xu, t, sum, ans;
    double a1 = a - 1.0, b1 = b - 1.0, mu = a / (a + b);
    double lnmu = log(mu), lnmuc = log(1. - mu);
    t = sqrt(a * b / (SQR(a + b) * (a + b + 1.0)));
    if (x > a / (a + b)) {
        if (x >= 1.0)
            return 1.0;
        xu = fmin(1., fmax(mu + 10. * t, x + 5.0 * t));
    } else {
        if (x <= 0.0)
            return 0.0;
        xu = fmax(0., fmin(mu - 10. * t, x - 5.0 * t));
    }
    sum = 0;
    for (j = 0; j < 18; j++) {
        t = x + (xu - x) * gauy[j];
        sum += gauw[j] * exp(a1 * (log(t) - lnmu) + b1 * (log(1 - t) - lnmuc));
    }
    ans = sum * (xu - x) *
          exp(a1 * lnmu - lgammaf(a) + b1 * lnmuc - lgammaf(b) + lgammaf(a + b));
    return ans > 0.0 ? 1.0 - ans : -ans;
}

double invbetai(double p, double a, double b) {
    const double EPS = 1.e-8;
    double pp, t, u, err, x, al, h, w, afac, a1 = a - 1., b1 = b - 1.;
    int j;
    if (p <= 0.)
        return 0.;
    else if (p >= 1.)
        return 1.;
    else if (a >= 1. && b >= 1.) {
        pp = (p < 0.5) ? p : 1. - p;
        t = sqrt(-2. * log(pp));
        x = (2.30753 + t * 0.27061) / (1. + t * (0.99229 + t * 0.04481)) - t;
        if (p < 0.5)
            x = -x;
        al = (SQR(x) - 3.) / 6.;
        h = 2. / (1. / (2. * a - 1.) + 1. / (2. * b - 1.));
        w = (x * sqrt(al + h) / h) - (1. / (2. * b - 1) - 1. / (2. * a - 1.)) *
                                         (al + 5. / 6. - 2. / (3. * h));
        x = a / (a + b * exp(2. * w));
    } else {
        double lna = log(a / (a + b)), lnb = log(b / (a + b));
        t = exp(a * lna) / a;
        u = exp(b * lnb) / b;
        w = t + u;
        if (p < t / w)
            x = pow(a * w * p, 1. / a);
        else
            x = 1. - pow(b * w * (1. - p), 1. / b);
    }
    afac = -lgammaf(a) - lgammaf(b) + lgammaf(a + b);
    //#pragma unroll
    for (j = 0; j < 10; j++) {
        if (x == 0. || x == 1.)
            return x;
        err = betai(a, b, x) - p;
        t = exp(a1 * log(x) + b1 * log(1. - x) + afac);
        u = err / t;
        x -= (t = u / (1. - 0.5 * fmin(1., u * (a1 / x - b1 / (1. - x)))));
        if (x <= 0.)
            x = 0.5 * (x + t);
        if (x >= 1.)
            x = 0.5 * (x + t + 1.);
        if (fabs(t) < EPS * x && j > 0)
            break;
    }
    return x;
}

// const double Beta::EPS = std::numeric_limits<double>::epsilon();
// const double Beta::FPMIN = std::numeric_limits<double>::min()/EPS;

/*
int main(){
        double x = invbetai(0.68547, 8., 10.);
        printf("%8f\n",x);
        return 0;
}
*/

#endif // PAEAN_INCGAMMABETA_H

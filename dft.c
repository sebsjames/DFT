// Implementation of DFT & iDFT
#include <stdio.h>
#include <math.h>

//define PI.
static const float PI = M_PI; // It is safer to explicitly define the type of the constant (instead of using #define)

// Some sample from a periodic function sin(0.1x/2pi) + cos(0.3x/2pi)
#define N_SAMPLES 18
// In samples, I have made all your constants single precision constants (adding the trailing f)
float samples[N_SAMPLES] = {1.047440989f, 0.8968022467f, 0.6104249648f, 0.2787682579f, 0.0f, -0.1420395219f, -0.09668181641f, 0.1420395219f, 0.5336978409f, 1.44167884f, 1.760073511f, 1.878694865f, 1.760073511f, 1.414213562f, 0.8968022467f, 0.2975560347f, -0.2787682579f, -0.7345720591f};

//Frequency constant warren. Integer multiples of this frequency are used to calculate F(k).
static const float w = 2.0f * PI / N_SAMPLES;

//Define a structure to use for complex numbers
typedef struct { float re, im; } complex;

//spectrum buffer in which results will be stored (same size as samples[] for inverse DFT)
complex spectrum[N_SAMPLES];

//evaluates the prevalence of an integer multiple of the fundamental frequency in some sample data.
//ran N_SAMPLES times for the entire DFT vector F(k).
complex DFT(float k)
{
    complex result = { 0.0f, 0.0f };

    //summation to find F(k) for some freqeuncy w*k.
    for (int sample_n = 0; sample_n < N_SAMPLES; sample_n++)
    {
        float angle = -w * k * (float)sample_n;
        result.re += samples[sample_n] * cos(angle);
        result.im += samples[sample_n] * sin(angle);
    }

    //divide by the number of samples.
    result.re /= (float)N_SAMPLES;
    result.im /= (float)N_SAMPLES;

    return result; // return strength (ray) and phase (angle) as a complex number.
}

//Inverse DFT
float iDFT(int sample_n) // rebuild a sample using the spectrum[] F(k).
{
    float result = 0.0f;

    //Inverse of DFT as outlined in README.md
    for (int freq_n = 0; freq_n < N_SAMPLES; freq_n++)
    {
        float angle = w * (float)freq_n * (float)sample_n;
        result += spectrum[freq_n].re * cos(angle) - spectrum[freq_n].im * sin(angle);
        result += spectrum[freq_n].re * sin(angle) + spectrum[freq_n].im * cos(angle);
    }

    //returns x_n and is executed N_SAMPLES times to reconstruct the full sample.
    return result;
}

int main()
{
    printf("test DFT F(k) ...\n");
    for (int i = 0; i < N_SAMPLES; i++)
    {
        //compute the frequency spectrum of F(n) from n 0->N_SAMPLES. Then output both componants.
        spectrum[i] = DFT((float)i);
        printf("[%f + i×%f]\n", spectrum[i].re, spectrum[i].im);
    }

    printf("|F(k)| ...\n");
    for (int k = 0; k < N_SAMPLES; k++)
    {
        float result = 0.0f;
        result = sqrt(spectrum[k].re * spectrum[k].re + spectrum[k].im * spectrum[k].im);
        printf("%f\n", result);
    }

    printf("test iDFT f(x) ...\n");
    for (int j = 0; j < N_SAMPLES; j++)
    {
        printf("[%f]\n", iDFT(j));
    }
    // The results of iDFT should be the same as samples[s] [so can you write a test to prove this?]

    return 0;
}

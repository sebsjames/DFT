// Implementation of DFT & iDFT
#include <stdio.h>
#include <math.h>

#define PI M_PI

// Some sample from a periodic function sin(0.1x/2pi) + cos(0.3x/2pi)
#define N_SAMPLES 18
float samples[N_SAMPLES] = {1.047440989, 0.8968022467, 0.6104249648, 0.2787682579, 0, -0.1420395219, -0.09668181641, 0.1420395219, 0.5336978409,1.44167884,1.760073511,1.878694865, 1.760073511, 1.414213562, 0.8968022467, 0.2975560347, -0.2787682579,-0.7345720591};

//Frequency constant warren. Integer multiples of this frequency is used in the samle doe
static float w = -2.0 * PI / N_SAMPLES;

typedef struct { float re, im; } complex; // structure to use for complex numbers

// spectrum buffer in which resultults will be stored (same size as samples[] for inverse DFT)
complex spectrum[N_SAMPLES];

// evaluates the previelince of an integer multiple of the fundamental freqeuncy in some sample data.
// ran N_SAMPLES times for the entire DFT vector.
complex DFT(float freq)
{
    complex result = { 0.0, 0.0 };

    //summation to find F(n) for some multiple of the fundamental freqeuncy n.
    for (int sample_n = 0; sample_n < N_SAMPLES; sample_n++)
    {
        float angle = w * freq * (float)sample_n;
        result.re += samples[sample_n] * cos(angle);
        result.im += samples[sample_n] * sin(angle);
    }

    //divide by the number of samples.
    result.re /= (float)N_SAMPLES;
    result.im /= (float)N_SAMPLES;

    return result; // return strength (ray) and phase (angle)
}

float iDFT(int sample_n) // rebuild a sample using the spectrum[]
{
    float result = 0.0;

    for (int freq_n = 0; freq_n < N_SAMPLES; freq_n++)
    {
        //result = sum(cosx - isinx + sinx + icosx) from samples 1->n

        float angle = 2.0 * PI * (float)freq_n * (float)sample_n / N_SAMPLES;
        result += spectrum[freq_n].re * cos(angle) - spectrum[freq_n].im * sin(angle);
        result += spectrum[freq_n].re * sin(angle) + spectrum[freq_n].im * cos(angle);
    }

    return result;
}

int main()
{
    printf("test DFT F(n) ...\n");
    for (int i = 0; i < N_SAMPLES; i++)
    {
        spectrum[i] = DFT((float)i);
        printf("[%f + i×%f]\n", spectrum[i].re, spectrum[i].im);
    }

    printf("|F(n)| ...\n");
    for (int k = 0; k < N_SAMPLES; k++)
    {
        float result = 0;
        result = sqrt(spectrum[k].re * spectrum[k].re + spectrum[k].im * spectrum[k].im);
        printf("%f\n", result);

    }

    printf("test iDFT ...\n");
    for (int j = 0; j < N_SAMPLES; j++)
    {
        printf("[%f]\n", iDFT(j));
    }
    // The resultults of iDFT should be the same as samples[s]

    return 0;
}

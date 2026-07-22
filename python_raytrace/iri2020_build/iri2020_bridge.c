#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(__GNUC__)
#define EXPORT __attribute__((visibility("default")))
#else
#define EXPORT
#endif

#define MAX_NUM_HEIGHTS 1000
#define OUTF_ROWS 20
#define OUTF_ROWS_KEEP 15
#define OARR_LEN 100

void iri2020_calc_(int *jf, int *jmag, float *glat, float *glon, int *year,
                   int *mmdd, float *dhour, float *heibeg, float *heiend,
                   float *heistp, float *outf, float *oarr);

static void set_error(char *error_msg, int error_len, const char *message) {
    if (error_msg == NULL || error_len <= 0) {
        return;
    }
    snprintf(error_msg, (size_t) error_len, "%s", message);
}

static void initialize_switches(int *jf) {
    int ii;
    for (ii = 0; ii < 3; ii++) jf[ii] = 1;
    for (ii = 3; ii < 6; ii++) jf[ii] = 0;
    for (ii = 6; ii < 20; ii++) jf[ii] = 1;
    for (ii = 20; ii < 23; ii++) jf[ii] = 0;
    for (ii = 23; ii < 27; ii++) jf[ii] = 1;
    for (ii = 27; ii < 30; ii++) jf[ii] = 0;
    for (ii = 30; ii < 32; ii++) jf[ii] = 1;
    jf[32] = 0;
    jf[33] = 1;
    jf[34] = 0;
    jf[35] = 1;
    jf[36] = 1;
    jf[37] = 1;
    jf[38] = 0;
    jf[39] = 0;
    jf[40] = 1;
    jf[41] = 1;
    jf[42] = 1;
    jf[43] = 1;
    jf[44] = 1;
    jf[45] = 1;
    jf[46] = 0;
    jf[47] = 1;
    jf[48] = 1;
    jf[49] = 1;

    jf[20] = 1;
    jf[21] = 0;
    jf[27] = 1;
    jf[33] = 0;
}

EXPORT int python_raytrace_iri2020_profile(
    double glat_deg,
    double glon_deg,
    double r12,
    int year,
    int month,
    int day,
    int hour,
    int minute,
    double heibeg_km,
    double heistp_km,
    int num_hts,
    int d_model_fpt,
    float *outf15,
    float *oarr_out,
    char *error_msg,
    int error_len
) {
    int jf[50];
    int jmag = 0;
    int mmdd;
    int ii;
    float glat;
    float glon;
    float heibeg;
    float heistp;
    float heiend;
    float dhour;
    float *outf_raw;
    float oarr_local[OARR_LEN];
    float ig12;
    float f107;

    if (num_hts < 2 || num_hts > MAX_NUM_HEIGHTS) {
        set_error(error_msg, error_len, "num_hts must be in the range [2, 1000]");
        return 1;
    }
    if (r12 <= 0.0 || r12 > 200.0) {
        set_error(error_msg, error_len, "R12 must be in the range (0, 200]");
        return 1;
    }
    if (heistp_km <= 0.0) {
        set_error(error_msg, error_len, "heistp_km must be positive");
        return 1;
    }

    initialize_switches(jf);
    jf[23] = d_model_fpt ? 0 : 1;

    jf[16] = 0;
    jf[24] = 0;
    jf[26] = 0;
    jf[31] = 0;
    jf[25] = 0;

    memset(oarr_local, 0, sizeof(oarr_local));
    f107 = (float) (63.75 + r12 * (0.728 + r12 * 0.00089));
    ig12 = (float) (-12.349154 + r12 * (1.4683266 - r12 * 2.67690893e-03));
    oarr_local[32] = (float) r12;
    oarr_local[38] = ig12;
    oarr_local[40] = f107;
    oarr_local[45] = f107;

    glat = (float) glat_deg;
    glon = (float) glon_deg;
    mmdd = month * 100 + day;
    dhour = (float) hour + ((float) minute / 60.0f) + 25.0f;
    heibeg = (float) heibeg_km;
    heistp = (float) heistp_km;
    heiend = heibeg + (float) (num_hts - 1) * heistp;

    outf_raw = (float *) calloc((size_t) (OUTF_ROWS * MAX_NUM_HEIGHTS), sizeof(float));
    if (outf_raw == NULL) {
        set_error(error_msg, error_len, "failed to allocate IRI2020 output buffer");
        return 1;
    }

    iri2020_calc_(jf, &jmag, &glat, &glon, &year, &mmdd, &dhour, &heibeg, &heiend, &heistp, outf_raw, oarr_local);

    for (ii = 0; ii < num_hts; ii++) {
        int jj;
        for (jj = 0; jj < OUTF_ROWS_KEEP; jj++) {
            outf15[ii * OUTF_ROWS_KEEP + jj] = outf_raw[ii * OUTF_ROWS + jj];
        }
    }
    for (ii = 0; ii < OARR_LEN; ii++) {
        oarr_out[ii] = oarr_local[ii];
    }

    free(outf_raw);
    return 0;
}

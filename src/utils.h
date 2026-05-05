static double normal_pdf_log(double x, double mean, double stddev) {
    double coefficient = 1.0 / (stddev * sqrt(2.0 * M_PI));
    double exponent = -(x - mean) * (x - mean) / (2 * stddev * stddev);
    return log(coefficient * exp(exponent));
}
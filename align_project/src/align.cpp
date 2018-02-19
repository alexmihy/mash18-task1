#include "align.h"
#include <string>

using std::string;
using std::cout;
using std::endl;
using std::tuple;
using std::make_tuple;
using std::tie;
using std::max;
using std::min;

#define MSEdRow 15
#define MSEdCol 15

long double calculateMSE(
    Image img1, Image img2, 
    int row1, int col1,
    int row2, int col2,
    int n_rows, int n_cols) 
{
    long long sum = 0;
    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            int r1, g1, b1;
            int r2, g2, b2;
            int c1, c2;

            tie(r1, g1, b1) = img1(row1 + i, col1 + j);
            tie(r2, g2, b2) = img2(row2 + i, col2 + j);

            c1 = (0.2126f * r1 + 0.7152f * g1 + 0.0722f * b1);
            c2 = (0.2126f * r2 + 0.7152f * g2 + 0.0722f * b2);

            sum += (c1 - c2) * (c1 - c2);
        }
    }

    long double result = sum / (1.0 * n_rows * n_cols);

    return result;
}

// TODO: functions align two images returns result

Image align(Image srcImage, bool isPostprocessing, std::string postprocessingType, double fraction, bool isMirror, 
            bool isInterp, bool isSubpixel, double subScale)
{
    int n_rows  = srcImage.n_rows;
    int n_cols = srcImage.n_cols;
    int out_n_rows = n_rows / 3;
    int out_n_cols = n_cols;

    Image outImage(out_n_rows, out_n_cols);

    long double MSE = -1.0;
    int MSE_row1 = 0, MSE_col1 = 0, MSE_row2 = 0, MSE_col2 = 0, MSE_n_rows = 0, MSE_n_cols = 0;
    for (int dRow = -MSEdRow; dRow <= MSEdRow; dRow++) {
        for (int dCol = -MSEdCol; dCol <= MSEdCol; dCol++) {
            int c_row1 = max(0, dRow); 
            int c_col1 = max(0, dCol);
            int c_row2 = min(0, abs(dRow)) + out_n_rows;
            int c_col2 = min(0, abs(dCol));
            int c_n_rows = out_n_rows - abs(dRow);
            int c_n_cols = out_n_cols - abs(dCol);


            long double cMSE = calculateMSE(srcImage, srcImage, c_row1, c_col1, c_row2, c_col2, c_n_rows, c_n_cols);

            if (MSE < 0 || cMSE < MSE) {
                MSE = cMSE;
                MSE_row1 = c_row1;
                MSE_col1 = c_col1;
                MSE_row2 = c_row2;
                MSE_col2 = c_col2;
                MSE_n_rows = c_n_rows;
                MSE_n_cols = c_n_cols;
            }
        }
    }

    printf("%d, %d, %d, %d, %d, %d\n", MSE_row1, MSE_col1, MSE_row2, MSE_col2, MSE_n_rows, MSE_n_cols);

    for (int i = 0; i < MSE_n_rows; i++) {
        for (int j = 0; j < MSE_n_cols; j++) {
            uint r, g, b;

            tie(r, r, r) = srcImage(i + MSE_row1, j + MSE_col1);
            tie(g, g, g) = srcImage(i + MSE_row2, j + MSE_col2);

            b = 0;

            outImage(i, j) = make_tuple(r, g, b);
        }
    }

    MSE = -1.0;
    for (int dRow = -MSEdRow; dRow <= MSEdRow; dRow++) {
        for (int dCol = -MSEdCol; dCol <= MSEdCol; dCol++) {
            int c_row1 = max(0, dRow); 
            int c_col1 = max(0, dCol);
            int c_row2 = min(0, abs(dRow)) + 2 * out_n_rows;
            int c_col2 = min(0, abs(dCol));
            int c_n_rows = out_n_rows - abs(dRow);
            int c_n_cols = out_n_cols - abs(dCol);


            long double cMSE = calculateMSE(outImage, srcImage, c_row1, c_col1, c_row2, c_col2, c_n_rows, c_n_cols);

            if (MSE < 0 || cMSE < MSE) {
                MSE = cMSE;
                MSE_row1 = c_row1;
                MSE_col1 = c_col1;
                MSE_row2 = c_row2;
                MSE_col2 = c_col2;
                MSE_n_rows = c_n_rows;
                MSE_n_cols = c_n_cols;
            }
        }
    }

    printf("%d, %d, %d, %d, %d, %d\n", MSE_row1, MSE_col1, MSE_row2, MSE_col2, MSE_n_rows, MSE_n_cols);

    for (int i = 0; i < MSE_n_rows; i++) {
        for (int j = 0; j < MSE_n_cols; j++) {
            uint r, g, b;

            tie(r, g, b) = outImage(i + MSE_row1, j + MSE_col1);
            tie(b, b, b) = srcImage(i + MSE_row2, j + MSE_col2);

            outImage(i, j) = make_tuple(r, g, b);
        }
    }

#if 0
    for (int i = 0; i < out_n_rows; i++) {
        for (int j = 0; j < out_n_cols; j++) {
            uint r, g, b;

            tie(r, r, r) = srcImage(i, j);
            tie(g, g, g) = srcImage(i + out_n_rows, j);
            tie(b, b, b) = srcImage(i + 2 * out_n_rows, j);

            outImage(i, j) = make_tuple(r, g, b);

        }
    }
#endif

    return outImage;
}

Image sobel_x(Image src_image) {
    Matrix<double> kernel = {{-1, 0, 1},
                             {-2, 0, 2},
                             {-1, 0, 1}};
    return custom(src_image, kernel);
}

Image sobel_y(Image src_image) {
    Matrix<double> kernel = {{ 1,  2,  1},
                             { 0,  0,  0},
                             {-1, -2, -1}};
    return custom(src_image, kernel);
}

Image unsharp(Image src_image) {
    return src_image;
}

Image gray_world(Image src_image) {
    return src_image;
}

Image resize(Image src_image, double scale) {
    return src_image;
}

Image custom(Image src_image, Matrix<double> kernel) {
    // Function custom is useful for making concrete linear filtrations
    // like gaussian or sobel. So, we assume that you implement custom
    // and then implement other filtrations using this function.
    // sobel_x and sobel_y are given as an example.
    return src_image;
}

Image autocontrast(Image src_image, double fraction) {
    return src_image;
}

Image gaussian(Image src_image, double sigma, int radius)  {
    return src_image;
}

Image gaussian_separable(Image src_image, double sigma, int radius) {
    return src_image;
}

Image median(Image src_image, int radius) {
    return src_image;
}

Image median_linear(Image src_image, int radius) {
    return src_image;
}

Image median_const(Image src_image, int radius) {
    return src_image;
}

Image canny(Image src_image, int threshold1, int threshold2) {
    return src_image;
}

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

#define mdRow 15
#define mdCol 15

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

long double calculateCC(
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

            sum += c1 * c2;
        }
    }

    long double result = sum;

    return result;
}

void align_two_image(
    Image &img1, Image &img2, 
    int mask1, int mask2, 
    Image &outImage, bool useMSE) 
{
    if (img1.n_rows != img2.n_rows || img1.n_cols != img2.n_cols)
        return;

    int n_rows  = img1.n_rows;
    int n_cols = img1.n_cols;

    long double metric = -1.0;
    int m_row1 = 0, m_col1 = 0, m_row2 = 0, m_col2 = 0, m_n_rows = 0, m_n_cols = 0;
    for (int dRow = -mdRow; dRow <= mdRow; dRow++) {
        for (int dCol = -mdCol; dCol <= mdCol; dCol++) {
            int c_row1 = max(0, dRow); 
            int c_col1 = max(0, dCol);
            int c_row2 = min(0, abs(dRow));
            int c_col2 = min(0, abs(dCol));
            int c_n_rows = n_rows - abs(dRow);
            int c_n_cols = n_cols - abs(dCol);

            long double c_metric;

            if (useMSE) {
                c_metric = calculateMSE(
                    img1, img2, 
                    c_row1, c_col1, c_row2, c_col2, c_n_rows, c_n_cols);
            } else {
                c_metric = calculateCC(
                    img1, img2, 
                    c_row1, c_col1, c_row2, c_col2, c_n_rows, c_n_cols);
            }

            if (metric < 0 || 
                (useMSE  && c_metric < metric) ||
                (!useMSE && c_metric > metric)) {
                metric = c_metric;
                m_row1 = c_row1;
                m_col1 = c_col1;
                m_row2 = c_row2;
                m_col2 = c_col2;
                m_n_rows = c_n_rows;
                m_n_cols = c_n_cols;
            }
        }
    }

    for (int i = 0; i < m_n_rows; i++) {
        for (int j = 0; j < m_n_cols; j++) {
            uint r1, g1, b1, r2, g2, b2;

            tie(r1, g1, b1) = img1(i + m_row1, j + m_col1);
            tie(r2, g2, b2) = img2(i + m_row2, j + m_col2);

            uint r, g, b;
            r = g = b = 0;

            if (mask1 & 4)
                r = r1;
            if (mask1 & 2)
                g = g1;
            if (mask1 & 1)
                b = b1;

            if (mask2 & 4)
                r = r2;
            if (mask2 & 2)
                g = g2;
            if (mask2 & 1)
                b = b2;

            outImage(i, j) = make_tuple(r, g, b);
        }
    }
}

Image align(Image srcImage, bool isPostprocessing, std::string postprocessingType, double fraction, bool isMirror, 
            bool isInterp, bool isSubpixel, double subScale)
{
    int n_rows  = srcImage.n_rows;
    int n_cols = srcImage.n_cols;
    int out_n_rows = n_rows / 3;
    int out_n_cols = n_cols;

    Image rImage(out_n_rows, out_n_cols);
    Image gImage(out_n_rows, out_n_cols);
    Image bImage(out_n_rows, out_n_cols);

    for (int k = 0; k < 3; k++) {
        for (int i = 0; i < out_n_rows; i++) {
            for (int j = 0; j < out_n_cols; j++) {
                uint r, g, b, br;

                tie(r, g, b) = srcImage(k * out_n_rows + i, j);
                br = (0.2126f * r + 0.7152f * g + 0.0722f * b);

                switch (k) {
                case 0:
                    rImage(i, j) = make_tuple(br, br, br);
                    break;
                case 1:
                    gImage(i, j) = make_tuple(br, br, br);
                    break;
                case 2:
                    bImage(i, j) = make_tuple(br, br, br);
                    break;
                default:
                    break;
                }
            }
        }
    }

    Image rgImage(out_n_rows, out_n_cols);
    Image outImage(out_n_rows, out_n_cols);

    align_two_image(rImage, gImage, 4, 2, rgImage, true); // 100, 010 => 110
    align_two_image(rgImage, bImage, 6, 1, outImage, true); // 110, 001 => 111

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

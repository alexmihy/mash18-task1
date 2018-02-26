#include "align.h"
#include <string>

#define _USE_MATH_DEFINES
#include <cmath> 

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

            c1 = (0.2126 * r1 + 0.7152 * g1 + 0.0722 * b1);
            c2 = (0.2126 * r2 + 0.7152 * g2 + 0.0722 * b2);

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

            c1 = (0.2126 * r1 + 0.7152 * g1 + 0.0722 * b1);
            c2 = (0.2126 * r2 + 0.7152 * g2 + 0.0722 * b2);

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
    #if 0
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
                br = (0.2126 * r + 0.7152 * g + 0.0722 * b);

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
    #endif

    //Image outImage = gaussian(srcImage, 1.4, 2);
    Image outImage = canny(srcImage, 1, 2);

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
    int n_rows = src_image.n_rows;
    int n_cols = src_image.n_cols;
    int rad = kernel.n_rows / 2;

    Image extImage(n_rows + 2 * rad, n_cols + 2 * rad);

    for (int i = 0; i < n_rows; i++)
        for (int j = 0; j < n_cols; j++)
            extImage(i + rad, j + rad) = src_image(i, j);

    for (int k = 0; k < rad; k++) {
        for (int i = 0; i < n_rows; i++) {
            extImage(i + rad, k) = src_image(i, rad - k);
            extImage(i + rad, n_cols + rad + k) = src_image(i, (n_cols - 1) - k);
        }

        for (int j = 0; j < n_cols; j++) {
            extImage(k, j + rad) = src_image(rad - k, j);
            extImage(n_rows + rad + k, j + rad) = src_image((n_rows - 1) - k, j);
        }
    }

    for (int i = 0; i < rad; i++) {
        for (int j = 0; j < rad; j++) {
            extImage(i, j) = src_image(rad - i, rad - j);
            extImage(i + n_rows + rad, j) = src_image((n_rows - 1) - i, rad - j);
            extImage(i, j + n_cols + rad) = src_image(rad - i, (n_cols - 1) - j);
            extImage(i + n_rows + rad, j + n_cols + rad) = src_image((n_rows - 1) - i, (n_cols - 1) - j);
        }
    }

    Image outImage(n_rows, n_cols);

    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            int R, G, B;
            double r_sum, g_sum, b_sum;

            r_sum = g_sum = b_sum = 0;
            for (int k1 = -rad; k1 <= rad; k1++) {
                for (int k2 = -rad; k2 <= rad; k2++) {
                    uint r, g, b;

                    tie(r, g, b) = extImage(i + rad + k1, j + rad + k2);
                    r_sum += r * kernel(k1 + rad, k2 + rad);
                    g_sum += g * kernel(k1 + rad, k2 + rad);
                    b_sum += b * kernel(k1 + rad, k2 + rad);
                }   
            }

            R = r_sum;
            G = g_sum; 
            B = b_sum;

            outImage(i, j) = make_tuple(
                max(min(R, 255), 0), 
                max(min(G, 255), 0), 
                max(min(B, 255), 0));
        }
    }

    return outImage;
}

Image autocontrast(Image src_image, double fraction) {
    return src_image;
}

Image gaussian(Image src_image, double sigma, int radius)  {
    int kernel_size = 2 * radius + 1;
    Matrix<double> kernel(kernel_size, kernel_size);

    long double sum = 0.0;
    for (int i = -radius; i <= radius; i++) {
        for (int j = -radius; j <= radius; j++) {
            double value;

            value = 1.0 / (2 * M_PI * sigma * sigma * exp((i * i + j * j) / (2 * sigma * sigma)));
            sum = sum + value;

            kernel(radius + i, radius + j) = value;
        }
    }

    double koef = 1 / sum;
    for (int i = 0; i < kernel_size; ++i)
        for (int j = 0; j < kernel_size; ++j)
            kernel(i, j) *= koef;

    return custom(src_image, kernel);
}

/* TODO: separable guassian */

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
    int n_rows  = src_image.n_rows;
    int n_cols = src_image.n_cols;

    Image bluredImage;
    Image Kx, Ky, K(n_rows, n_cols);

    Image grayImg(n_rows, n_cols);

    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            uint r, g, b, br;

            tie(r, g, b) = src_image(i, j);
            br = (0.2126 * r + 0.7152 * g + 0.0722 * b);

            grayImg(i, j) = make_tuple(br, br, br);
        }
    }

    bluredImage = gaussian(grayImg, 1.4, 2);

    Kx = sobel_x(bluredImage);
    Ky = sobel_y(bluredImage);

    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            uint b1, b2, bret;

            tie(b1, b1, b1) = Kx(i, j);
            tie(b2, b2, b2) = Ky(i, j);

            bret = sqrt(b1 * b1 + b2 * b2);

            K(i, j) = make_tuple(bret, bret, bret);
        }
    }



    return K;
}

#include "pch.h"
#include"mySparseLibrary.h"
#include "kingghidorah_S.h"
#include "kingghidorah_C.h"

void main()
{
    int n = 1;
    int s = 10;
    bool parallel = true;
    int M = 100 * n;
    int N = 200 * n;
    auto mat = gcnew kingghidorah::mySparse(M, N);
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j += s) {
            mat->adddat(i, j, 1);
        }
    }
    mat->merge();
    mat->clearcoeff();
    if (parallel)
        mat->ofAtA(mat);
    else
        mat->_ofAtA(mat);
    std::cout<<mat->_at(0, 0) << std::endl;
    std::cout<<mat->num_elem(0)<<std::endl;
    mat = gcnew kingghidorah::mySparse(M, N);
    mat->begin_construct();
    for (int ii = 0; ii < 10; ii++)
    {
        auto _mat = gcnew kingghidorah::mySparse(M, N);
        for (int i = 0; i < M / 10; i++)
        {
            for (int j = 0; j < N; j += s) {
                _mat->adddat(i + M * ii / 10, j, 1);

                _mat->addcoeff(1.0);
            }
        }
        mat->addmat(_mat);
    }
    mat->resize(mat->rows(), mat->cols());
    mat->merge();
    std::cout<<mat->rows() << std::endl;
    std::cout << mat->cols() << std::endl;
    std::cout << mat->num_elem(3) << std::endl;
    mat->freezecoeff();
    if (parallel)
        mat->ofAtA(mat);
    else
        mat->_ofAtA(mat);
    std::cout<<mat->_at(0, 0) << std::endl;
    std::cin.get();
}
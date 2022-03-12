#include"mySparseLibrary.h"
void KingOfMonsters::helper::contract(System::Collections::Generic::List<workspace^>^ _mats, denseMatrix^ U, denseMatrix^ V, denseMatrix^ W, array<sparseMatrix^>^ _mats1, array<sparseMatrix^>^ _mats2, array<sparseMatrix^>^ _mats3, int _C)
{
	int _mt = omp_get_max_threads();
	int m = _mats->Count;
	int n = U->get().cols();
	//array<denseMatrix^>^ mm = gcnew array<mySparse^>(m);
	std::vector<Eigen::SparseMatrix<double>> mm01(m); //multiplied with U
	std::vector<Eigen::SparseMatrix<double>> mm02(m); //multiplied with V
	std::vector<Eigen::SparseMatrix<double>> mm03(m); //multiplied with V
	array<sparseMatrix^>^ mm11 = gcnew array<sparseMatrix^>(n + 1); //mutiplied with U and W
	array<sparseMatrix^>^ mm12 = gcnew array<sparseMatrix^>(n + 1); //mutiplied with V and W
	array<sparseMatrix^>^ mm13 = gcnew array<sparseMatrix^>(n + 2); //mutiplied with V and W

	int count = 0;
	//denseMatrix^ Uo = gcnew denseMatrix();
	//Uo->get() = U->get() * (U->get().transpose() * U->get()).inverse();
	//denseMatrix^ Vo = gcnew denseMatrix();
	//Vo->get() = U->get() * (U->get().transpose() * U->get()).inverse();
#pragma omp parallel
	{
#pragma omp single
		_mt = omp_get_num_threads();
	}
#pragma omp parallel for
	for (int ii = 0; ii < _mt; ii++)
	{
		int S = ii * m / _mt;
		int E = (ii + 1) * m / _mt;
		for (int i = S; i < E; i++)
		{
			auto M = _mats[i];
			auto newM1 = M->tosparse(n, true, false);
			auto newM2 = M->tosparse(n, false, false);
			auto newM3 = M->tosparse(n, false, true);
			mm01[i] = newM1;
			mm02[i] = newM2;
			mm03[i] = newM3;
#pragma omp critical
			{
				count++;
				Console::WriteLine(count.ToString() + "/" + m.ToString());
			}
		}
	}

	//array<KingOfMonsters::denseMatrix^>::Resize(_mats1, m);
	//array<KingOfMonsters::denseMatrix^>::Resize(_mats2, m);
	//array<KingOfMonsters::denseMatrix^>::Resize(_mats3, m);

	
	for (int i = 0; i < m; i++)
	{
		_mats1[i] = gcnew sparseMatrix(n, n);
		_mats2[i] = gcnew sparseMatrix(n, n);
		_mats3[i] = gcnew sparseMatrix(n, n);
		_mats1[i]->get() = mm01[i];
		_mats2[i]->get() = mm02[i];
		_mats3[i]->get() = mm03[i];
	}
	/*auto _W = W->get();
	count = 0;
#pragma omp parallel for
	for (int ii = 0; ii < _mt; ii++)
	{
		int S = ii * n / _mt;
		int E = (ii + 1) * n / _mt;
		for (int i = S; i < E; i++)
		{
			//auto M = _mats[i];

			mm11[i] = gcnew denseMatrix(n, n);
			mm12[i] = gcnew denseMatrix(n, n);
			mm13[i] = gcnew denseMatrix(n, n);

			mm11[i]->get().setZero();
			mm12[i]->get().setZero();
			mm13[i]->get().setZero();
			for (int k = 0; k < m; k++)
			{
				if (mm01[k].cols() == n && mm01[k].rows() == n)
				{
					mm11[i]->get() += mm01[k] * _W(k, i);
					mm12[i]->get() += mm02[k] * _W(k, i);
					mm13[i]->get() += mm03[k] * _W(k, i);
				}
			}
#pragma omp critical
			{
				count++;
				Console::WriteLine(count.ToString() + "/" + n.ToString());
			}
		}
	}

	mm11[n] = gcnew denseMatrix(n, n);
	mm11[n]->get().setZero();
	mm13[n] = gcnew denseMatrix(n, n);
	mm13[n]->get().setZero();
	mm13[n + 1] = gcnew denseMatrix(n, n);
	mm13[n + 1]->get().setZero();

	for (int i = 0; i < n; i++)
	{
		for (int k = 0; k < m; k++)
		{
			if (mm01[k].rows() == 1 && mm01[k].cols() == n)
			{
				mm11[n]->get().row(i) += mm01[k] * _W(k, i);
				mm13[n]->get().row(i) += (mm03[k] * _W(k, i)).transpose();
			}
		}
	}

	mm12[n] = gcnew denseMatrix(n, n);
	mm12[n]->get().setZero();
	for (int i = 0; i < n; i++)
	{
		for (int k = 0; k < m; k++)
		{
			if (mm02[k].rows() == 1 && mm02[k].cols() == n)
			{
				mm12[n]->get().row(i) += mm02[k] * _W(k, i);
				mm13[n + 1]->get().row(i) += (mm03[k] * _W(k, i));
			}
		}
	}
	for (int i = 0; i < n; i++)
	{
		_mats1[i] = mm11[i];
		_mats2[i] = mm12[i];
		_mats3[i] = mm13[i];
	}
	_mats1[n] = mm11[n];
	_mats2[n] = mm12[n];
	_mats3[n] = mm13[n];////wu
	_mats3[n + 1] = mm13[n + 1];//wv
	*/

}
double KingOfMonsters::helper::VarPro(myDoubleArray^ phi, myDoubleArray^ zz, denseMatrix^ __U, denseMatrix^ __V, denseMatrix^ __W, array<sparseMatrix^>^ _mats1, array<sparseMatrix^>^ _mats2, array<sparseMatrix^>^ _mats3, myDoubleArray^ _r1, myDoubleArray^ _r2, double dt,int tt)
{
	std::vector<Eigen::SparseMatrix<double>*> __mats1;
	std::vector<Eigen::SparseMatrix<double>*> __mats2;
	std::vector<Eigen::SparseMatrix<double>*> __mats3;
	__mats1.clear();
	__mats2.clear();
	__mats3.clear();

	for (int i = 0; i < _mats1->Length; i++)
	{
		__mats1.push_back(&_mats1[i]->get());
	}
	for (int i = 0; i < _mats2->Length; i++)
	{
		__mats2.push_back(&_mats2[i]->get());
	}
	for (int i = 0; i < _mats3->Length; i++)
	{
		__mats3.push_back(&_mats3[i]->get());
	}
	double norm=KingOfMonsters::_helper::VarPro(&phi->_arr->__v, &zz->_arr->__v, &__U->get(), &__V->get(), &__W->get(), __mats1, __mats2, __mats3, &_r1->_arr->__v, &_r2->_arr->__v, dt, tt);
	System::Console::WriteLine("normGrad=" + norm.ToString());
	return norm;
}
double KingOfMonsters::helper::ALT(myDoubleArray^ phi, myDoubleArray^ zz, denseMatrix^ __U, denseMatrix^ __V, denseMatrix^ __W, array<sparseMatrix^>^ _mats1, array<sparseMatrix^>^ _mats2, array<sparseMatrix^>^ _mats3, myDoubleArray^ _r1, myDoubleArray^ _r2, double dt, int tt)
{
	std::vector<Eigen::SparseMatrix<double>*> __mats1;
	std::vector<Eigen::SparseMatrix<double>*> __mats2;
	std::vector<Eigen::SparseMatrix<double>*> __mats3;
	__mats1.clear();
	__mats2.clear();
	__mats3.clear();

	for (int i = 0; i < _mats1->Length; i++)
	{
		__mats1.push_back(&_mats1[i]->get());
	}
	for (int i = 0; i < _mats2->Length; i++)
	{
		__mats2.push_back(&_mats2[i]->get());
	}
	for (int i = 0; i < _mats3->Length; i++)
	{
		__mats3.push_back(&_mats3[i]->get());
	}
	double norm=KingOfMonsters::_helper::ALT(&phi->_arr->__v, &zz->_arr->__v, &__U->get(), &__V->get(), &__W->get(), __mats1, __mats2, __mats3, &_r1->_arr->__v, &_r2->_arr->__v, dt, tt);
	System::Console::WriteLine("normGrad="+norm.ToString());
	return norm;
}
void KingOfMonsters::helper::write(myDoubleArray^ phi0, myDoubleArray^ zz0,myDoubleArray^ phi, myDoubleArray^ zz, denseMatrix^ __U, denseMatrix^ __V, denseMatrix^ __W, array<sparseMatrix^>^ _mats1, array<sparseMatrix^>^ _mats2, array<sparseMatrix^>^ _mats3, myDoubleArray^ _r1, myDoubleArray^ _r2, double dt, int tt)
{
	System::Console::WriteLine(System::Environment::CurrentDirectory);

	std::vector<Eigen::SparseMatrix<double>*> __mats1;
	std::vector<Eigen::SparseMatrix<double>*> __mats2;
	std::vector<Eigen::SparseMatrix<double>*> __mats3;
	__mats1.clear();
	__mats2.clear();
	__mats3.clear();

	for (int i = 0; i < _mats1->Length; i++)
	{
		__mats1.push_back(&_mats1[i]->get());
	}
	for (int i = 0; i < _mats2->Length; i++)
	{
		__mats2.push_back(&_mats2[i]->get());
	}
	for (int i = 0; i < _mats3->Length; i++)
	{
		__mats3.push_back(&_mats3[i]->get());
	}
	KingOfMonsters::_helper::write(&phi0->_arr->__v, &zz0->_arr->__v ,&phi->_arr->__v, &zz->_arr->__v, &__U->get(), &__V->get(), &__W->get(), __mats1, __mats2, __mats3, &_r1->_arr->__v, &_r2->_arr->__v, dt, tt);


}
void KingOfMonsters::helper::computeKrylovSubspace(System::Collections::Generic::List<workspace^>^ _mats, denseMatrix^ _U, denseMatrix^ _V, denseMatrix^ _W,  int nU, int nV, int r, myPermutation^ mphi, myPermutation^ mZ, myDoubleArray^ phi, myDoubleArray^ zz, System::Collections::Generic::List<Tuple<int, int>^>^ bb1, System::Collections::Generic::List<Tuple<int, int>^>^bb2)
{
	int _C = phi->_arr->__v.size();
	if (r > _C)r = _C;
	int m = _mats->Count;
	Eigen::MatrixXd U(_C, r);
	Eigen::MatrixXd V(_C, r);
	Eigen::MatrixXd W(m, r);

	Eigen::VectorXd u0 = phi->_arr->__v;
	Eigen::VectorXd v0 = zz->_arr->__v;
	//u0.applyOnTheLeft(mphi->p->perm);
	//v0.applyOnTheLeft(mZ->p->perm);
	//System::Collections::Generic::List<Tuple<int, int>^>^ bb1 = gcnew System::Collections::Generic::List<Tuple<int, int>^>();
	//System::Collections::Generic::List<Tuple<int, int>^>^ bb2 = gcnew System::Collections::Generic::List<Tuple<int, int>^>();

	for (int i = 0; i < _C; i++)
	{
		if (mphi->p->perm.indices()(i) >= nU)
		{
			workspace^ work = gcnew workspace();
			work->_dat->push_back(_Triplet<double>(i, -1, 1));
			_mats->Add(work);
			bb1->Add(gcnew Tuple<int, int>(i,_mats->Count-1));
		}
	}
	for (int i = 0; i < _C; i++)
	{
		if (mZ->p->perm.indices()(i) >= nV)
		{
			workspace^ work = gcnew workspace();
			work->_dat->push_back(_Triplet<double>(-1, i, 1));
			_mats->Add(work);
			bb2->Add(gcnew Tuple<int, int>(i, _mats->Count - 1));
		}
	}


	/*Eigen::VectorXd w0(m);
	Eigen::VectorXd ui(_C);
	Eigen::VectorXd vi(_C);
	Eigen::VectorXd wi(m);

	ui.setRandom();
	vi.setRandom();
	wi.setRandom();
	ui = u0;
	vi = v0;
	//wi.setOnes();

	
	ui.normalize();
	vi.normalize();
	wi.normalize();


	U.setZero();
	V.setZero();
	W.setZero();
	U.col(0) = ui;
	V.col(0) = vi;
	W.col(0) = wi;
	U.col(0).normalize();
	V.col(0).normalize();
	W.col(0).normalize();

	u0 = U.col(0);
	v0 = V.col(0);
	w0 = W.col(0);

	Eigen::VectorXd u2(_C);
	Eigen::VectorXd v2(_C);
	Eigen::VectorXd w2(m);
	u2.setZero();
	v2.setZero();
	w2.setZero();
	int _mt = omp_get_max_threads();
	#pragma omp parallel
	{
	#pragma omp single
		_mt = omp_get_num_threads();
	}

	Console::WriteLine("task#:" + (0).ToString() + "/" + r.ToString());
	int pointer = 0;
	for (int j = 0; j < r - 1; j++)
	{
		//ui.topRows(nU) = U.col(j);
		//vi.topRows(nV) = V.col(j);
		//ui.bottomRows(_C - nU) = u0.bottomRows(_C - nU);
		//vi.bottomRows(_C - nV) = v0.bottomRows(_C - nV);
		ui = U.col(j);
		vi = V.col(j);
		wi = W.col(j);

		u2.setZero();
		v2.setZero();
		w2.setZero();

		u0 = U.col(pointer);
		v0 = V.col(pointer);
		w0 = W.col(pointer);
		//if (j > (pointer + 1) * 3)
		//{
		//	pointer++;
		//}

		Console::WriteLine("task#:" + (j + 1).ToString() + "/" + r.ToString());
		//int count = 0;
#pragma omp parallel for
		for (int ii = 0; ii < _mt; ii++)
		{
			int S = ii * m / _mt;
			int E = (ii + 1) * m / _mt;
			Eigen::VectorXd _u2(_C);
			Eigen::VectorXd _v2(_C);
			_u2.setZero();
			_v2.setZero();
			//double* ptr = &w2.coeffRef(S);
			//double* ptr2 = &wi.coeffRef(S);
			for (int i = S; i < E; i++)
			{

				auto M = _mats[i];

				w2.coeffRef(i) = M->mulboth(&u0, &vi,j);// (ui.transpose() * M)* vi;
				//ptr++;

				M->mulright(&v0, &_u2, wi.coeffRef(i), j);
				M->mulleft(&ui, &_v2, wi.coeffRef(i), j);// (ui.transpose() * M).transpose();

			}
#pragma omp critical
			{
				u2 += _u2;
				v2 += _v2;
			}
		}
		for (int kk = 0; kk < 10; kk++)
		{
			Eigen::VectorXd hu = U.leftCols(j + 1).transpose() * u2;
			Eigen::VectorXd hv = V.leftCols(j + 1).transpose() * v2;
			Eigen::VectorXd hw = W.leftCols(j + 1).transpose() * w2;
			if (hu.norm() < 0.0000000000001)
			{
				if (hv.norm() < 0.0000000000001)
				{
					if (hw.norm() < 0.0000000000001)
					{
						break;
					}
				}
			}
			u2 = u2 - U.leftCols(j + 1) * hu;
			v2 = v2 - V.leftCols(j + 1) * hv;
			w2 = w2 - W.leftCols(j + 1) * hw;

			u2.normalize();
			v2.normalize();
			w2.normalize();

		}
		U.col(j + 1) = u2;
		V.col(j + 1) = v2;
		W.col(j + 1) = w2;
	}*/
	_U->resize(_C, r);
	_V->resize(_C, r);
	_W->resize(m, r);
	_U->get().setZero();
	_V->get().setZero();
	_W->get().setZero();
	U.setIdentity();
	V.setIdentity();
	_U->get() = U;
	_V->get() = V;
	_W->get() = W;
}

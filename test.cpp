#include "four4angel.h"
#include "math.h"
#include<iostream>
using namespace std;

// ��� ����� ������������


four4angel::four4angel(void)
{
	for (int i = 0; i < NUM_POINTS;i++)
	{
		p[i].x = 0.0;
		p[i].y = 0.0;
		cout << "I don't have problems" << endl;
	}
}
void four4angel::calcVectors() // ôóíêöèÿ ñ÷èòàþùàÿ è âûâîäÿùàÿ âåêòîðà 
{
	p[NUM_POINTS].x = p[0].x;
	p[NUM_POINTS].y = p[0].y;
	for (int i = 0; i < NUM_POINTS;i++)
	{
		v[i].x = p[i + 1].x - p[i].x;
		v[i].y = p[i + 1].y - p[i].y;
		cout << "Êîîðäèíàòû " << i+1 << " âåêòîðà " << v[i].x << '\t' << v[i].y << endl;
	}
	d1.x = p[2].x - p[0].x;
	d1.y = p[2].y - p[0].y;
	d2.x = p[1].x - p[3].x;
	d2.y = p[1].y - p[3].y;
}

four4angel::four4angel(point* _p)
{
	for (int i = 0; i < NUM_POINTS;i++)
	{
		p[i].x = _p[i].x;
		p[i].y = _p[i].y;
	}
	calcVectors();
}
four4angel::four4angel(vect* _v)
{
	p[0].x = 0.0;
	p[0].y = 0.0;
	for (int i = 0; i < NUM_POINTS-1; i++)
	{
		v[i].x = _v[i].x;
		v[i].y = _v[i].y;
	}
	for (int i = 1; i < NUM_POINTS; i++)
	{
		p[i].x = p[i - 1].x + v[i - 1].x;
		p[i].y = p[i - 1].y + v[i - 1].y;
	}
	calcVectors();
}
four4angel::four4angel(double _s1, double _s2, double _s3, double _s4, double _u)
{
	_u = _u * Pi / 180; // ïåðåâîäèì â ðàäèàíû
	p[0].x = 0.0;
	p[0].y = 0.0;
	p[1].x = _s1;
	p[1].y = 0.0;
	p[3].x = _s4 * round(cos(_u));
	p[3].y = _s4 * round(sin(_u));
	double Diagonal_naprotiv_ugla = sqrt(pow(_s4, 2) + pow(_s1, 2) - 2 *_s1*_s4* cos(_u)); // ïî òåîðåìå êîñèíóñîâ
	double Ugol_4_1 = acos((pow(_s4, 2) - pow(_s1, 2) - pow(Diagonal_naprotiv_ugla, 2)) / ((-2) * Diagonal_naprotiv_ugla * _s1)); // óãîë ïîä äèàãîíàëüþ ñâåðõó
	double Ugol_4_2= acos((pow(_s3, 2) - pow(_s2, 2) - pow(Diagonal_naprotiv_ugla, 2)) / ((-2) * _s2 * Diagonal_naprotiv_ugla)); // óãîë íàä äèàãîíàëüþ ñâåðõó
	double Ugol_4 = Ugol_4_1 + Ugol_4_2; 
	double Diagonal_iz_ugla = sqrt(pow(_s1, 2) + pow(_s2, 2) - 2 * _s1 * _s2 * cos(Ugol_4)); // ïî òåîðåìå êîñèíóñîâ
	double Ugol_1_1 = asin((_s2 * sin(Ugol_4)) / (Diagonal_iz_ugla)); // ÷àñòü èçâåñòíîãî óãëà ïîä äèàãîíàëüþ, âûõîäÿùåé èç íåãî

	p[2].x = Diagonal_iz_ugla * cos(Ugol_1_1);
	p[2].y = Diagonal_iz_ugla * sin(Ugol_1_1);

	calcVectors();
}
void four4angel::calcVectors3d() // ôóíêöèÿ ñ÷èòàþùàÿ è âûâîäÿùàÿ âåêòîðà 
{
	p[NUM_POINTS].x = p[0].x;
	p[NUM_POINTS].y = p[0].y;
	for (int i = 0; i < NUM_POINTS;i++)
	{
		v[i].x = p[i + 1].x - p[i].x;
		v[i].y = p[i + 1].y - p[i].y;
		v[i].z = 0.0;
		v2[i].x = 0.0;
		v2[i].y = 0.0;
		v2[i].z = h;
		cout << "Êîîðäèíàòû " << i + 1 << " âåêòîðà " << v[i].x << '\t' << v[i].y << '\t' << v[i].z << endl;
	}
	for (int i = 0; i < NUM_POINTS;i++)
	{
		cout << "Êîîðäèíàòû " << i + 1 + 4 << " âåêòîðà " << v2[i].x << '\t' << v2[i].y << '\t' << v2[i].z << endl;
	}
	for (int i = 0; i < NUM_POINTS;i++)
	{
		cout << "Êîîðäèíàòû " << i + 1 + 8 << " âåêòîðà " << v[i].x << '\t' << v[i].y << '\t' << v[i].z << endl;
	}
	d1.x = p[2].x - p[0].x;
	d1.y = p[2].y - p[0].y;
	d2.x = p[1].x - p[3].x;
	d2.y = p[1].y - p[3].y;
}
four4angel::four4angel(point* _p, double _h)
{
	h = _h;
	for (int i = 0; i < NUM_POINTS;i++)
	{
		p[i].x = _p[i].x;
		p[i].y = _p[i].y;
		p[i].z = 0.0;
		p2[i].x = _p[i].x;
		p2[i].y = _p[i].y;
		p2[i].z = h;
	}
	calcVectors3d();
}
four4angel::four4angel(vect* _v, double _h)
{
	h = _h;
	p[0].x = 0.0;
	p[0].y = 0.0;
	p[0].z = 0.0;
	for (int i = 0; i < NUM_POINTS - 1; i++)
	{
		v[i].x = _v[i].x;
		v[i].y = _v[i].y;
		v[i].z = 0.0;
	}
	for (int i = 1; i < NUM_POINTS; i++)
	{
		p[i].x = p[i - 1].x + v[i - 1].x;
		p[i].y = p[i - 1].y + v[i - 1].y;
		p[i].z = 0.0;
		p2[i].x = p[i - 1].x + v[i - 1].x;
		p2[i].y = p[i - 1].y + v[i - 1].y;
		p2[i].z = h;
	}
	calcVectors3d();
}
four4angel::four4angel(double _s1, double _s2, double _s3, double _s4, double _u, double _h)
{
	h = _h;
	_u = _u * Pi / 180; // ïåðåâîäèì â ðàäèàíû
	p[0].x = p2[0].x = 0.0;
	p[0].y = p2[0].y = 0.0;
	p[0].z = 0.0;
	p2[0].z = h;
	p[1].x = p2[1].x = _s1;
	p[1].y = p2[1].y = 0.0;
	p[1].z = 0.0;
	p2[1].z = h;
	p[3].x = p2[3].x = _s4 * round(cos(_u));
	p[3].y = p2[3].y = _s4 * round(sin(_u));
	p[3].z = 0.0;
	p2[3].z = h;
	double Diagonal_naprotiv_ugla = sqrt(pow(_s4, 2) + pow(_s1, 2) - 2 * _s1 * _s4 * cos(_u)); // ïî òåîðåìå êîñèíóñîâ
	double Ugol_4_1 = acos((pow(_s4, 2) - pow(_s1, 2) - pow(Diagonal_naprotiv_ugla, 2)) / ((-2) * Diagonal_naprotiv_ugla * _s1)); // óãîë ïîä äèàãîíàëüþ ñâåðõó
	double Ugol_4_2 = acos((pow(_s3, 2) - pow(_s2, 2) - pow(Diagonal_naprotiv_ugla, 2)) / ((-2) * _s2 * Diagonal_naprotiv_ugla)); // óãîë íàä äèàãîíàëüþ ñâåðõó
	double Ugol_4 = Ugol_4_1 + Ugol_4_2;
	double Diagonal_iz_ugla = sqrt(pow(_s1, 2) + pow(_s2, 2) - 2 * _s1 * _s2 * cos(Ugol_4)); // ïî òåîðåìå êîñèíóñîâ
	double Ugol_1_1 = asin((_s2 * sin(Ugol_4)) / (Diagonal_iz_ugla)); // ÷àñòü èçâåñòíîãî óãëà ïîä äèàãîíàëüþ, âûõîäÿùåé èç íåãî

	p[2].x = p2[2].x = Diagonal_iz_ugla * cos(Ugol_1_1);
	p[2].y = p2[2].y = Diagonal_iz_ugla * sin(Ugol_1_1);
	p[2].z = 0.0;
	p2[2].z = h;

	calcVectors3d();
}
double scalVect(point v1, point v2) // ñêàëÿðíîå ïðîèçâåäåíèå âåêòîðîâ
{
	return v1.x * v2.x + v1.y * v2.y;
};
double ugolvect(point v1, point v2) // íàõîæäåíèå óãëà â ðàäèàíàõ ìåæäó äâóìÿ âåêòîðàìè
{
	return acos(scalVect(v1, v2) / sqrt(scalVect(v1, v1) * scalVect(v2, v2)));
};
void four4angel::PrintCord(int i) // âûâîä êîîðäèíàò òî÷åê
{
	if (i-1 >= 4||i-1<0)
		cout << "Òàêîé òî÷êè íåò"<< endl;
	else
		cout << "Êîîðäèíàòû " << i << " òî÷êè " << "  " << p[i-1].x << '\t' << p[i-1].y << endl;
}
void four4angel::PrintCord3d(int i) // âûâîä êîîðäèíàò òî÷åê
{
	if (i - 1 >= 8 || i - 1 < 0)
		cout << "Òàêîé òî÷êè íåò" << endl;
	else
		if (i <= 4)
			cout << "Êîîðäèíàòû " << i << " òî÷êè " << "  " << p[i - 1].x << '\t' << p[i - 1].y << '\t' << p[i - 1].z << endl;
		else
			cout << "Êîîðäèíàòû " << i << " òî÷êè " << "  " << p2[i - 5].x << '\t' << p2[i - 5].y << '\t' << p2[i - 5].z << endl;
}
void four4angel::PrintDlin(int i) // âûâîä äëèí ñòîðîí
{
	if (i - 1 >= 4 || i - 1 < 0)
		cout << "Òàêîé ñòîðîíû íåò";
	else
		cout << "Äëèíà "<< i << " ñòîðîíû ðàâíà "<< mv[i-1];
}
void four4angel::PrintDlin3d(int i) // âûâîä äëèí ñòîðîí
{
	if (i - 1 >= 12 || i - 1 < 0)
		cout << "Òàêîé ñòîðîíû íåò";
	else
		if(i<=4)
			cout << "Äëèíà " << i << " ñòîðîíû ðàâíà " << mv[i - 1];
	    if (i>4 && i<=8 )
		    cout << "Äëèíà " << i << " ñòîðîíû ðàâíà " << h;
		if (i>8)
			cout << "Äëèíà " << i << " ñòîðîíû ðàâíà " << mv[i - 9];
}
void four4angel::Diagonali() // âûâîä äëèí äèàãîíàëåé
{
	d11 = sqrt(pow((p[0].x - p[2].x), 2) + pow((p[0].y - p[2].y), 2));
	d22 = sqrt(pow((p[1].x - p[3].x), 2) + pow((p[1].y - p[3].y), 2));
	cout << d11 << '\t' << d22;
}
void four4angel::Diagonali3d() // âûâîä äëèí äèàãîíàëåé
{
	d11 = sqrt(pow((p[0].x - p[2].x), 2) + pow((p[0].y - p[2].y), 2));
	d22 = sqrt(pow((p[1].x - p[3].x), 2) + pow((p[1].y - p[3].y), 2));
	d3d1 = sqrt(pow(d11, 2) + pow(h, 2));
	d3d2 = sqrt(pow(d22, 2) + pow(h, 2));
	cout << "Äèàãîíàëè îñíîâàíèé ðàâíû " << d11 << '\t' << d22<< endl;
	cout << "Äëèíû äèàãîíàëåé ïðèçìû ðàâíû " << d3d1 << '\t' << d3d2;
}
double four4angel::Perimetr() // âû÷èñëåíèå ïåðèìåòðà
{
	for (int i = 0; i < NUM_POINTS;i++)
	{
		mv[i] = sqrt(v[i].x * v[i].x + v[i].y * v[i].y);
		Per += mv[i];
	}
	return Per;
}
double four4angel::Perimetr3d() // âû÷èñëåíèå ïåðèìåòðà
{
	for (int i = 0; i < NUM_POINTS;i++)
	{
		mv[i] = sqrt(v[i].x * v[i].x + v[i].y * v[i].y);
		Per += mv[i];
	}
	return Per*2 + 4*h;
}
double four4angel::Ploshad() // âû÷èñëåíèå ïëîùàäè
{
	S=d11*d22*sin(ugolvect(d1,d2))/2;
	return S;
}
double four4angel::Ploshad3d() // âû÷èñëåíèå ïëîùàäè ïîëíîé ïîâåðõíîñòè
{
	return d11 * d22 * sin(ugolvect(d1, d2)) + h * mv[0] + h * mv[1] + h * mv[2] + h * mv[3];
}
double four4angel::Volume()
{
	return (d11 * d22 * sin(ugolvect(d1, d2)) / 2) * h;
}
double four4angel::PrintUgol(int i) // âûâîä óãëà ïî íîìåðó âåðøèíû
{
	ugol[0] = Pi - ugolvect(v[0], v[3]);
	ugol[1] = Pi - ugolvect(v[0], v[1]);
	ugol[2] = Pi - ugolvect(v[1], v[2]);
	ugol[3] = Pi - ugolvect(v[2], v[3]);

	for (int i = 0; i < NUM_POINTS;i++)
	{
		if (ugol[i] > Pi / 2)
			ugol[i] = 2 * Pi - ugol[i];
	}
	if (i - 1 >= 4 || i - 1 < 0)
		cout << "Òàêîé âåðøèíû íåò";
	else
		cout << "Óãîë ïðè " << i << " âåðøèíå ðàâåí " << ugol[i - 1] << " ðàäèàí èëè æå ïðèáëèçèòåëüíî " << round(ugol[i - 1] * 180 / Pi) << " ãðàäóñîâ " << endl;
	return ugol[i - 1];
}
double four4angel::PrintUgol(int i, int a) // âûâîä óãëà ïî ñòîðîíàì ìåæäó êîòîðûìè îí ëåæèò
{
	double r = Pi - ugolvect(v[i - 1], v[a - 1]);
	if (r > Pi / 2)
		r = 2 * Pi - r;
	if (i - 1 >= 4 || i - 1 < 0 || -1 >= 4 || a - 1 < 0)
		cout << "Òàêîé ñòîðîíû íåò" << endl;
	else
		cout << "Óãîë ìåæäó " << i << " è " << a << " ñòîðîíàìè ðàâåí  " << r << " ðàäèàí èëè æå ïðèáëèçèòåëüíî " << round(r * 180 / Pi) << " ãðàäóñîâ " << endl;
	return r;
}
void four4angel::Vogn() // ïðîâåðêà íà âîãíóòîñòü/âûïóêëîñòü
{
	d1.x = p[2].x - p[0].x;
	d1.y = p[2].y - p[0].y;
	d2.x = p[1].x - p[3].x;
	d2.y = p[1].y - p[3].y;
	double v1 = d2.x * v[3].y - d2.y * v[3].x;
	double v2 = d2.x * (p[2].y - p[3].y) - d2.y * (p[2].x - p[3].x);
	double v3 = d1.x * (p[3].y - p[0].y) - d1.y * (p[3].x - p[0].x);
	double v4 = d1.x * v[0].y - d1.y * v[0].x;
	if (v1 * v2 < 0 && v4 * v3 < 0)
		cout << "×åòûð¸õóãîëüíèê âûïóêëûé";
	else
		cout << "×åòûð¸õóãîëüíèê âîãíóòûé";
}
void four4angel::Vogn3d() // ïðîâåðêà íà âîãíóòîñòü/âûïóêëîñòü
{
	d1.x = p[2].x - p[0].x;
	d1.y = p[2].y - p[0].y;
	d2.x = p[1].x - p[3].x;
	d2.y = p[1].y - p[3].y;
	double v1 = d2.x * v[3].y - d2.y * v[3].x;
	double v2 = d2.x * (p[2].y - p[3].y) - d2.y * (p[2].x - p[3].x);
	double v3 = d1.x * (p[3].y - p[0].y) - d1.y * (p[3].x - p[0].x);
	double v4 = d1.x * v[0].y - d1.y * v[0].x;
	if (v1 * v2 < 0 && v4 * v3 < 0)
		cout << "Îñíîâàíèå ïðèçìû âûïóêëîå";
	else
		cout << "Îñíîâàíèå ïðèçìû âîãíóòîå";
		cout << "Ýòî ÿ äîáàâèëà äëÿ ãèòà";
}
four4angel::~four4angel(void) // äåñòðóêòîð
{

}

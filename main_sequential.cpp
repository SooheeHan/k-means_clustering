#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

//#include "windows.h"
//#include "omp.h"

// variables
int m_nRow, m_nColumn, m_nBand, m_nClass, m_nMaxIter;
double m_nChangedRatio;
unsigned short **p_Image;
unsigned char *p_Class;
double **p_Centers;

int *p_nMembers;

timeval t_start0, t_end0;
timeval t_start, t_end;

int main()
{
	// Initializaing
	
	/*m_nRow=7841;
	m_nColumn=7691;
	m_nBand=7;
	m_nClass=10;
	m_nMaxIter=100;
	m_nChangedRatio=0.015;

	char inPath[100], outPath[100], statPath[100];
	sprintf(inPath, "/mirror/data/landsat8.bin");
	sprintf(outPath, "/mirror/proc/landsat8_classified_squential.bin");
	sprintf(statPath, "/mirror/proc/landsat8_classified_squential.txt");*/

	/*m_nRow=30978;
	m_nColumn=40980;
	m_nBand=8;
	m_nClass=10;
	m_nMaxIter=100;
	m_nChangedRatio=0.020;

	char inPath[100], outPath[100], statPath[100];
	sprintf(inPath, "/mirror/data/S2A_20160518T094844.bin");
	sprintf(outPath, "/mirror/proc/S2A_20160518T094844_classified_squential.bin");
	sprintf(statPath, "/mirror/proc/S2A_20160518T094844_classified_squential.txt");*/

	/*m_nRow=29080;
	m_nColumn=30320;
	m_nBand=4;
	m_nClass=10;
	m_nMaxIter=20;
	m_nChangedRatio=0.015;
	char * inPath = "D:/Data/multi_sharpened.bin";
	char * outPath = "D:/Data/multi_sharpened_classified_squential.bin";
	char * statPath = "D:/Data/multi_sharpened_classified_squential.txt";*/

	m_nRow=7051;
	m_nColumn=1761;
	m_nBand=242;
	m_nClass=10;
	m_nMaxIter=20;
	m_nChangedRatio=0.015;
	char inPath[100], outPath[100], statPath[100];
	sprintf(inPath, "/mirror/data/eo.bin");
	sprintf(outPath, "/mirror/proc/eo_classified_squential.bin");
	sprintf(statPath, "/mirror/proc/eo_classified_squential.txt");

	FILE *pf;
	pf = fopen(statPath,"wt");

	p_Image = new unsigned short*[m_nBand];	//p_Image[m_nBand]�� Ŭ����
	
	// Reading image
	//t_start0 = gettimeofday(&t_start, 0);
	gettimeofday(&t_start0, 0);
	gettimeofday(&t_start, 0);

	FILE *p_file = fopen(inPath,"rb");

	for(int i=0; i<m_nBand; i++)
	{
		p_Image[i] = new unsigned short[m_nRow*m_nColumn];
		fread(p_Image[i], sizeof(unsigned short), m_nRow*m_nColumn, p_file);
	}

	fclose(p_file);

	p_Class = new unsigned char[m_nRow*m_nColumn];
	for (int i=0; i<m_nRow*m_nColumn; i++)
		p_Class[i]=254;
		
	gettimeofday(&t_end, 0);
	printf ("Image is read in %.2f seconds\n", float(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_usec-t_start.tv_usec)/1000000.0);
	fprintf (pf, "Image is read in %.2f seconds\n", float(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_usec-t_start.tv_usec)/1000000.0);


	// �ʱ� Ŭ���� �߽� ����
	gettimeofday(&t_start, 0);

	int nEffectivePixel = m_nRow*m_nColumn;
	p_Centers = new double *[m_nClass];
	p_nMembers = new int[m_nClass];

	double **p_Centers_temp = new double*[m_nClass];
	for(int i=0; i<m_nClass; i++)
	{
		p_Centers[i] = new double[m_nBand];
		p_Centers_temp[i] = new double[m_nBand];
		p_nMembers[i] = 0;
	}

	for (int i=0; i<m_nBand; i++)
	{
		double min, max;
		min = 1000000;
		max = -1000000;

		for (int j=0; j<m_nRow*m_nColumn; j++)
		{
			if (p_Image[i][j]==0)	//NULL value skip
			{
				if (p_Class[j] == 254)
				{
					p_Class[j] = 255;
					nEffectivePixel --; //Null pixel�� ��ü ī��Ʈ���� ����
				}
				continue;
			}

			if (p_Image[i][j]<min)
				min=p_Image[i][j];
			if (p_Image[i][j]>max)
				max=p_Image[i][j];
		}

		double interval = (max-min)/double(m_nClass);

		for (int j=0; j<m_nClass; j++)
		{
			p_Centers[j][i] = min+j*interval+interval/2;
			fprintf(pf, "band %d class %d min %.4lf max %.4lf center %.4lf\n", i, j, min, max, p_Centers[j][i]);
		}
	}
	gettimeofday(&t_end, 0);
	printf ("Class is initialized in %.2f seconds\n", float(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_usec-t_start.tv_usec)/1000000.0);
	fprintf (pf, "Class is initialized in %.2f seconds\n", float(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_usec-t_start.tv_usec)/1000000.0);

	for (int iter=0; iter<m_nMaxIter; iter++)
	{
		for(int i=0; i<m_nClass; i++)
		{
			for (int j=0; j<m_nBand; j++)
				p_Centers_temp[i][j] = 0;
			p_nMembers[i] = 0;
		}

		// �з� ����
		gettimeofday(&t_start, 0);
		int nChanged = 0;
		int nearest;

		for (int i=0; i<m_nRow*m_nColumn; i++)
		{
			if (p_Class[i] == 255)
				continue;

			double dist = 10000000000;
			for (int j=0; j<m_nClass; j++)
			{
				if (p_Centers[j][0] == -1)
					continue;

				double cur_dist = 0;
				for (int k=0; k<m_nBand; k++)
					cur_dist += pow(p_Image[k][i]-p_Centers[j][k],2);

				cur_dist = sqrt(cur_dist);
				if (cur_dist < dist)
				{
					dist = cur_dist;
					nearest = j;
				}
			}

			if (p_Class[i] != nearest)
			{
				nChanged ++;
				p_Class[i] = nearest;
			}

			int nClass = p_Class[i];
			for (int j=0; j<m_nBand; j++)
				p_Centers_temp[nClass][j] += p_Image[j][i];

			p_nMembers[nClass]++;
		}

		// Ŭ���� �߽� �缳��
		for (int i=0; i<m_nClass; i++)
			for (int j=0; j<m_nBand; j++)
			{
				if (p_nMembers[i] != 0)
					p_Centers[i][j] = p_Centers_temp[i][j]/p_nMembers[i];
				else
					p_Centers[i][j] = -1;
			}

		gettimeofday(&t_end, 0);
		double nChangedRatio = double(nChanged)/double(nEffectivePixel);
		printf ("%d th iteration: Classified in %.2f seconds and %.4lf changed\n", iter, float(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_usec-t_start.tv_usec)/1000000.0, nChangedRatio);
		fprintf (pf, "%d th iteration: Classified in %.2f seconds and %.4lf changed\n", iter, float(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_usec-t_start.tv_usec)/1000000.0, nChangedRatio);

		// stat ���� -----------------------------------------------------------------------------
		/*
				for (int i=0; i<m_nClass; i++)
							for (int j=0; j<m_nBand; j++)
								fprintf(pf, "class %d band %d center %.2lf member %d \n", i, j, p_Centers[i][j], p_nMembers[i]);*/
				

		if (nChangedRatio<m_nChangedRatio)
			break;
	}

	// Writing result
	gettimeofday(&t_start, 0);

	p_file = fopen(outPath,"wb");
	fwrite(p_Class, sizeof(unsigned char), m_nRow*m_nColumn, p_file);
	fclose(p_file);

	gettimeofday(&t_end, 0);
	printf ("Classification result is written in %.2f seconds\n", float(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_usec-t_start.tv_usec)/1000000.0);
	fprintf (pf, "Classification result is written in %.2f seconds\n", float(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_usec-t_start.tv_usec)/1000000.0);

	// Releasing array
	for(int i=0; i<m_nBand; i++)
		delete p_Image[i];
	delete p_Image;

	for(int i=0; i<m_nClass; i++)
	{
		delete p_Centers[i];
		delete p_Centers_temp[i];
	}
	delete p_Centers;
	delete p_Centers_temp;

	delete p_nMembers;

	//t_end0 = GetTickCount();
	gettimeofday(&t_end0, 0);

	printf ("The whole process is done in %.2f seconds\n", float(t_end0.tv_sec-t_start0.tv_sec)+(t_end0.tv_usec-t_start0.tv_usec)/1000000.0);
	fprintf (pf, "The whole process is done in %.2f seconds\n", float(t_end0.tv_sec-t_start0.tv_sec)+(t_end0.tv_usec-t_start0.tv_usec)/1000000.0);

	fclose(pf);

	return 0;
	//system("pause");
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "omp.h"
#include <sys/time.h>

// variables
int m_nRow, m_nColumn, m_nBand, m_nClass, m_nMaxIter;
double m_nChangedRatio;
unsigned short **p_Image;
unsigned char *p_Class;
double **p_Centers;
double ***p_Centers_thread;

int *p_nMembers;
int **p_nMembers_thread;

timeval t_start0, t_end0;
timeval t_start, t_end;

int main()
{
	int m_nThread;
	m_nThread = omp_get_num_procs();

	// Initializaing
	
	// for the sample data, Landsat8.bin
	m_nRow=7841;
	m_nColumn=7691;
	m_nBand=7;
	m_nClass=10;
	m_nMaxIter=100;
	m_nChangedRatio=0.015;
	
	char inPath[100], outPath[100], statPath[100];
	sprintf(inPath, "/mirror/data/landsat8.bin");
	sprintf(outPath, "/mirror/proc/landsat8_classified_OpenMP.bin");
	sprintf(statPath, "/mirror/proc/landsat8_classified_OpenMP.txt");

	// for the sample data, sentinel2a.bin (not included to the repository because of too large volume)
	/*m_nRow=30978;
	m_nColumn=40980;
	m_nBand=8;
	m_nClass=10;
	m_nMaxIter=20;
	m_nChangedRatio=0.015;

	char inPath[100], outPath[100], statPath[100];
	sprintf(inPath, "/mirror/data/S2A_20160518T094844.bin");
	sprintf(outPath, "/mirror/proc/S2A_20160518T094844_classified_OpenMP.bin");
	sprintf(statPath, "/mirror/proc/S2A_20160518T094844_classified_OpenMP.txt");*/
	
	FILE *pf;
	pf = fopen(statPath,"wt");

	p_Image = new unsigned short*[m_nBand];
	
	// Reading image
	//t_start0 = t_start = GetTickCount();
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

	 
	// declaring global variables for class centers 
	gettimeofday(&t_start, 0);
	p_Centers = new double *[m_nClass];
	p_nMembers = new int[m_nClass];

	for(int i=0; i<m_nClass; i++)
	{
		p_Centers[i] = new double[m_nBand];
		p_nMembers[i] = 0;
	}

	// declaring local variables for class centers  
	p_Centers_thread = new double **[m_nThread];
	p_nMembers_thread = new int*[m_nThread];
	for(int i=0; i<m_nThread; i++)
	{
		p_Centers_thread[i] = new double *[m_nClass];
		p_nMembers_thread[i] = new int[m_nClass];

		for(int j=0; j<m_nClass; j++)
		{
			p_Centers_thread[i][j] = new double[m_nBand];
			p_nMembers_thread[i][j] = 0;
		}
	}

	// calculating initial class centers 
	int nEffectivePixel = m_nRow*m_nColumn;
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
					nEffectivePixel --; //don't count Null pixels
				p_Class[j] = 255;
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
	fflush(stderr);

	fprintf (pf, "Class is initialized in %.2f seconds\n", float(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_usec-t_start.tv_usec)/1000000.0);


	// classification loop
	for (int iter=0; iter<m_nMaxIter; iter++)
	{
		gettimeofday(&t_start, 0);
		int i, j, k, nCurThread, nClass;
		int nChanged = 0;
		double dist; 
		double cur_dist; 
		int nearest;

		// initializing local class centers
#pragma omp parallel private(nCurThread, i, j)
		{
			nCurThread = omp_get_thread_num();
			for (i=0; i<m_nClass; i++)
			{
				for (j=0; j<m_nBand; j++)
					p_Centers_thread[nCurThread][i][j] = 0;

				p_nMembers_thread[nCurThread][i] = 0;
			}
		}

		// classification
#pragma omp parallel for private(i, j, k, dist, cur_dist, nearest, nClass, nCurThread) reduction(+:nChanged)
		for (i=0; i<m_nRow*m_nColumn; i++)
		{
			if (p_Class[i] == 255)
				continue;

			dist = 10000000000;
			for (j=0; j<m_nClass; j++)
			{
				if (p_Centers[j][0] == -1)
					continue;

				cur_dist = 0;
				for (k=0; k<m_nBand; k++)				
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

			nCurThread = omp_get_thread_num();
			nClass = p_Class[i];

			for (j=0; j<m_nBand; j++)
				p_Centers_thread[nCurThread][nClass][j] += p_Image[j][i];

			p_nMembers_thread[nCurThread][nClass]++;
		}

		// recalculating class centers
		for (int i=0; i<m_nClass; i++)
			for (int j=0; j<m_nBand; j++)
			{
				p_Centers[i][j] = 0;
				p_nMembers[i] = 0;

				for (nCurThread = 0; nCurThread<m_nThread; nCurThread++)
				{
					p_Centers[i][j] += p_Centers_thread[nCurThread][i][j];
					p_nMembers[i] += p_nMembers_thread[nCurThread][i];
				}

				if (p_nMembers[i] != 0)
					p_Centers[i][j] = p_Centers[i][j]/p_nMembers[i];
				else
					p_Centers[i][j] = -1;
			}

		gettimeofday(&t_end, 0);
		double nChangedRatio = double(nChanged)/double(nEffectivePixel);
		printf ("%d th iteration: Classified in %.2f seconds and %.4lf changed\n", iter, float(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_usec-t_start.tv_usec)/1000000.0, nChangedRatio);
		fprintf (pf, "%d th iteration: Classified in %.2f seconds and %.4lf changed\n", iter, float(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_usec-t_start.tv_usec)/1000000.0, nChangedRatio);

		// writing stats		
		for (int i=0; i<m_nClass; i++)
			for (int j=0; j<m_nBand; j++)
				fprintf(pf, "class %d band %d center %.2lf member %d \n", i, j, p_Centers[i][j], p_nMembers[i]);
		
		if (nChangedRatio<m_nChangedRatio)
			break;
	}
	
	gettimeofday(&t_start, 0);

	p_file = fopen(outPath,"w+b");
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
		delete p_Centers[i];
	delete p_Centers;

	delete p_nMembers;

	for (int nCurThread = 0; nCurThread<m_nThread; nCurThread++)
	{
		for(int i=0; i<m_nClass; i++)
			delete p_Centers_thread[nCurThread][i];

		delete p_nMembers_thread[nCurThread];
	}
	delete p_Centers_thread;
	delete p_nMembers_thread;	

	gettimeofday(&t_end0, 0);
	printf ("The whole process is done in %.2f seconds\n", float(t_end0.tv_sec-t_start0.tv_sec)+(t_end0.tv_usec-t_start0.tv_usec)/1000000.0);
	fprintf (pf, "The whole process is done in %.2f seconds\n", float(t_end0.tv_sec-t_start0.tv_sec)+(t_end0.tv_usec-t_start0.tv_usec)/1000000.0);

	fclose(pf);
	return 0;
}

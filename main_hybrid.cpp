#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "omp.h"
#include <sys/time.h>
//#include <unistd.h>
//#include <ctime>

// variables
int m_nRow, m_nColumn, m_nBand, m_nClass, m_nMaxIter;
double m_nChangedRatio;
int m_nPixelsInNode, m_nPixelsPerNode;
int m_nEffectivePixel;
int m_bIteration;
int m_nSelect;
int m_nBuffer;

unsigned short **p_Image;
unsigned char *p_Class;
double **p_Centers;
int *p_nMembers;

timeval t_start0, t_end0;
timeval t_start, t_end;
timeval temp_start, temp_end;

using namespace std;

int main(int argc, char* argv[])
{
	// MPI variables

	int m_nProcessor;
	int m_nID;
	MPI_Status status;

	// MPI implementation

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&m_nProcessor);
	MPI_Comm_rank(MPI_COMM_WORLD,&m_nID);
	char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);


	fprintf(stderr, "Using processor %s, rank %d out of %d processors\n", processor_name, m_nID, m_nProcessor);

	int m_nThread;
	m_nThread = omp_get_num_procs();

	// Initializing

	// for the sample data, Landsat8.bin
	m_nRow=7841;
	m_nColumn=7691;
	m_nBand=7;
	m_nClass=10;
	m_nMaxIter=100;
	m_nChangedRatio=0.015;

	char inPath[100], outPath[100], localPath[100], statPath[100], transPath[100];
	sprintf(inPath, "/mirror/data/landsat8.bin"); // change it to a correct path in a master node
	sprintf(outPath, "/mirror/proc/landsat8_classified_MPI_%d.bin", m_nProcessor-1); // change it to a correct path in a master node
	sprintf(localPath, "/home/cluster/landsat8.bin"); // change it to a correct path in each slave node
	sprintf(statPath, "/mirror/proc/landsat8_classified_MPI_%d.txt", m_nProcessor-1); // change it to a correct path in a master node
	sprintf(transPath, "/mirror/proc/landsat8_transfer.txt"); // change it to a correct path in a master node
	
	// for the sample data, sentinel2a.bin (not included to the repository because of too large volume)
	/*m_nRow=30978;
	m_nColumn=40980;
	m_nBand=8;
	m_nClass=10;
	m_nMaxIter=20;
	m_nChangedRatio=0.015;

	char inPath[100], outPath[100], localPath[100], statPath[100], transPath[100];
	sprintf(inPath, "/mirror/data/S2A_20160518T094844.bin");
	sprintf(outPath, "/mirror/proc/S2A_20160518T094844_classified_MPI_%d.bin", m_nProcessor-1);
	sprintf(localPath, "/home/cluster/S2A_20160518T094844.bin");
	sprintf(statPath, "/mirror/proc/S2A_20160518T094844_classified_MPI_%d.txt", m_nProcessor-1);
	sprintf(transPath, "/mirror/proc/S2A_20160518T094844_transfer.txt");*/

	FILE *pf;	

	// 0. declaration of variables	
	m_nPixelsPerNode = m_nPixelsInNode = ceil(double(m_nRow*m_nColumn)/double(m_nProcessor-1));
	if (m_nID == m_nProcessor-1)
		m_nPixelsInNode = m_nRow*m_nColumn-(m_nProcessor-2)*m_nPixelsPerNode;

	m_nEffectivePixel = m_nRow*m_nColumn;
	m_bIteration = 1;
	m_nBuffer = 8000000;	// buffer size of image transfer from master to nodes using MPI

	if (m_nID == 0) // master node
	{
		// selecting task
		fprintf(stderr, "Select case 1.create local file, 2.use local file, 3.use remote file : ");
		fflush(stderr);
		scanf("%d", &m_nSelect);
	}

	gettimeofday(&t_start0, 0);

	MPI_Bcast(&m_nSelect, 1, MPI_INT, 0, MPI_COMM_WORLD);

	unsigned short *pBuffer;
	
	FILE *p_file;
	int nIter, nRemaining;

	switch(m_nSelect)
	{
		case 1:	// transferring image from master to nodes, and saving it locally in each node
			if (m_nID == 0)
				pf = fopen(transPath,"wt");

			gettimeofday(&t_start, 0);
			if (m_nID == 0)
			{
				fprintf(stderr, "Head : reading and transferring image\n");
				fflush(stderr);
				fprintf(pf, "Head : reading and transferring image\n");
				p_file = fopen(inPath,"rb");
			}
			else
			{
				p_file = fopen(localPath,"wb");
			}
							
			pBuffer = new unsigned short[m_nBuffer];
			nIter = int(double(m_nRow*m_nColumn)/double(m_nBuffer));
			nRemaining = m_nRow*m_nColumn - m_nBuffer*nIter;

			for(int i=0; i<m_nBand; i++)
			{
				for (int j=0; j<nIter; j++)
				{
					if (m_nID == 0)
						fread(pBuffer, sizeof(unsigned short), m_nBuffer, p_file);

					MPI_Bcast(pBuffer, m_nBuffer, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);

					if (m_nID != 0)
						fwrite(pBuffer, sizeof(unsigned short), m_nBuffer, p_file);
				}

				
				if (m_nID == 0)
						fread(pBuffer, sizeof(unsigned short), nRemaining, p_file);

				MPI_Bcast(pBuffer, nRemaining, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);

				if (m_nID != 0)
					fwrite(pBuffer, sizeof(unsigned short), nRemaining, p_file);
			}
			delete pBuffer;
			fclose(p_file);

			gettimeofday(&t_end, 0);
			if (m_nID == 0)
			{
				fprintf(stderr, "Head : image is transfered and local file created in %.2f seconds\n", float(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_usec-t_start.tv_usec)/1000000.0);
				fflush(stderr);
				fprintf(pf, "Head : image is transfered and local file created in %.2f seconds\n", float(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_usec-t_start.tv_usec)/1000000.0);
			}

			MPI_Finalize();
			return 0;



		case 2:	// processing using the locally saved file in each node
			if (m_nID == 0)
				pf = fopen(statPath,"wt");

			MPI_Barrier(MPI_COMM_WORLD);
			gettimeofday(&t_start, 0);
			if (m_nID != 0)
			{				
				// 1. loading image
				p_file = fopen(localPath, "rb");

				p_Image = new unsigned short*[m_nBand];
				for(int i=0; i<m_nBand; i++)
				{
					p_Image[i] = new unsigned short[m_nPixelsInNode];
					int64_t pos = (int64_t(i)*int64_t(m_nRow*m_nColumn) + int64_t(m_nID-1)*int64_t(m_nPixelsPerNode))*int64_t(sizeof(unsigned short));
					
					fseeko64(p_file, pos, SEEK_SET);	//64bit file position
					fread(p_Image[i], sizeof(unsigned short), m_nPixelsInNode, p_file);
				}

				fclose(p_file);
				gettimeofday(&t_end, 0);
				fprintf(stderr, "Node %d : local file is read in %.2f seconds\n", m_nID, float(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_usec-t_start.tv_usec)/1000000.0);
				fflush(stderr);

				p_Class = new unsigned char[m_nPixelsInNode];
				for (int i=0; i<m_nPixelsInNode; i++)
					p_Class[i]=254;
			}
			
			MPI_Barrier(MPI_COMM_WORLD);
			if (m_nID == 0)
			{
				gettimeofday(&t_end, 0);
				fprintf(stderr, "Head : local file is read in %.2f seconds\n", float(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_usec-t_start.tv_usec)/1000000.0);
				fflush(stderr);
				fprintf(pf, "Head : local file is read in %.2f seconds\n", float(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_usec-t_start.tv_usec)/1000000.0);
			}
			
			break;

		case 3:	// transferring pieces of image from master to nodes, and processing directly without saving it
			if (m_nID == 0)
				pf = fopen(statPath,"wt");

			gettimeofday(&t_start, 0);
			if (m_nID == 0) // master
			{
				// 1. distributing image from master to nodes
				fprintf(stderr, "Head : reading and transfering image\n");
				fflush(stderr);
				fprintf(pf, "Head : reading and transfering image\n");
	
				p_file = fopen(inPath,"rb");

				unsigned short * pBuffer = new unsigned short[m_nPixelsPerNode];
				for(int i=0; i<m_nBand; i++)
				{
					for (int j=1; j<=m_nProcessor-2; j++)
					{
						fread(pBuffer, sizeof(unsigned short), m_nPixelsPerNode, p_file);
						MPI_Send(pBuffer, m_nPixelsPerNode, MPI_UNSIGNED_SHORT, j, 0, MPI_COMM_WORLD);
					}

					int nRemaining = m_nRow*m_nColumn-(m_nProcessor-2)*m_nPixelsPerNode;
					fread(pBuffer, sizeof(unsigned short), nRemaining, p_file);
					MPI_Send(pBuffer, nRemaining, MPI_UNSIGNED_SHORT, m_nProcessor-1, 0, MPI_COMM_WORLD);
				}
				delete pBuffer;
				fclose(p_file);

				gettimeofday(&t_end, 0);
				fprintf(stderr, "Image is read and transfered in %.2f seconds\n", float(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_usec-t_start.tv_usec)/1000000.0);
				fflush(stderr);
				fprintf(pf, "Image is read and transfered in %.2f seconds\n", float(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_usec-t_start.tv_usec)/1000000.0);
			}
			else // nodes
			{				
				p_Image = new unsigned short*[m_nBand];
				for(int i=0; i<m_nBand; i++)
				{
					p_Image[i] = new unsigned short[m_nPixelsInNode];
					MPI_Recv(p_Image[i], m_nPixelsInNode, MPI_UNSIGNED_SHORT, 0, 0, MPI_COMM_WORLD, &status);	
				}

				fprintf(stderr, "Node%d : received %d pixels\n", m_nID, m_nPixelsInNode);
				fflush(stderr);

				p_Class = new unsigned char[m_nPixelsInNode];
				for (int i=0; i<m_nPixelsInNode; i++)
					p_Class[i]=254;
			}
			break;

		default:
			MPI_Finalize();
			return 0;
	}



	// 2. initializing class centers
	gettimeofday(&t_start, 0);

	p_Centers = new double *[m_nClass];
	for(int i=0; i<m_nClass; i++)
		p_Centers[i] = new double[m_nBand];
		
	p_nMembers = new int[m_nClass];

	int nNullPixel = 0;
	int nLocalNullPixel = 0;
	double nMax, nMin, nLocalMax, nLocalMin;

	for (int i=0; i<m_nBand; i++)
	{
		nLocalMin = 1000000;
		nLocalMax = -1000000;

		if (m_nID != 0)
		{
			for (int j=0; j<m_nPixelsInNode; j++)
			{
				if (p_Image[i][j]==0)	//NULL value skip
				{
					if (p_Class[j] == 254)
						nLocalNullPixel++;	//don't count Null pixes
					p_Class[j] = 255;
					continue;
				}
				if (nLocalMin > p_Image[i][j]) nLocalMin=p_Image[i][j];
				if (nLocalMax < p_Image[i][j]) nLocalMax=p_Image[i][j];
			}
		}

		MPI_Reduce(&nLocalMin, &nMin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
		MPI_Reduce(&nLocalMax, &nMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

		if (m_nID == 0)
		{
			double interval = (nMax-nMin)/double(m_nClass);

			for (int j=0; j<m_nClass; j++)
			{
				p_Centers[j][i] = nMin+double(j)*interval+interval/2;
				fprintf(pf, "band %d class %d min %.4lf max %.4lf center %.4lf\n", i, j, nMin, nMax, p_Centers[j][i]);
			}

		}
	}

	MPI_Reduce(&nLocalNullPixel, &nNullPixel, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	m_nEffectivePixel -= nNullPixel;

	// tranferring class centers
	for (int i=0; i<m_nClass; i++)
		MPI_Bcast(p_Centers[i], m_nBand, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	gettimeofday(&t_end, 0);
	if (m_nID == 0)
	{
		fprintf(stderr, "Head : class is initialized in %.2f seconds\n", float(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_usec-t_start.tv_usec)/1000000.0);
		fflush(stderr);
		fprintf(pf, "Head : class is initialized in %.2f seconds\n", float(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_usec-t_start.tv_usec)/1000000.0);
	}

	// a variable for master
	int nTotalChanged;

	// declaring global variables for class centers 
	double ** p_LocalCenters = new double *[m_nClass];
	for(int i=0; i<m_nClass; i++)
		p_LocalCenters[i] = new double[m_nBand];
	int *p_nLocalMembers = new int[m_nClass];
	
	// declaration of local variables for class centers 
	double ***p_Centers_thread;
	int **p_nMembers_thread;
	int nLocalChanged;

	if (m_nID != 0) // nodes
	{
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
	}

	// 3. classificaiton loop
	for (int iter=0; iter<m_nMaxIter; iter++)
	{
		gettimeofday(&t_start, 0); // redo
		nTotalChanged = nLocalChanged = 0;

		for (int i=0; i<m_nClass; i++)
		{
			for (int j=0; j<m_nBand; j++)
				p_LocalCenters[i][j] = 0;

			p_nLocalMembers[i] = 0;
		}

		if (m_nID != 0)
		{
			int nCurThread;
			int i, j, k;
			double dist, cur_dist; 
			int nearest, nClass;

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
			gettimeofday(&temp_start, 0);
#pragma omp parallel for private(i, j, k, dist, cur_dist, nearest, nClass, nCurThread) reduction(+:nLocalChanged)
			for (i=0; i<m_nPixelsInNode; i++)
			{
				if (p_Class[i] == 255)
					continue;

				dist = 10000000000; // max distance not to be reached
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
					nLocalChanged ++;
					p_Class[i] = nearest;
				}

				nCurThread = omp_get_thread_num();
				nClass = p_Class[i];

				for (j=0; j<m_nBand; j++)
					p_Centers_thread[nCurThread][nClass][j] += p_Image[j][i];

				p_nMembers_thread[nCurThread][nClass]++;
			}
			gettimeofday(&temp_end, 0);
			fprintf(stderr, "Node %d, running time %.2lf, nPixel %d, nClass %d, nBand %d\n", 
				m_nID, float(temp_end.tv_sec-temp_start.tv_sec)+(temp_end.tv_usec-temp_start.tv_usec)/1000000.0, m_nPixelsInNode, m_nClass, m_nBand); fflush(stderr); //temporary

			// recalculating class centers in nodes
			for (int i=0; i<m_nClass; i++)
				for (nCurThread = 0; nCurThread<m_nThread; nCurThread++)				
				{
					for (int j=0; j<m_nBand; j++)
						p_LocalCenters[i][j] += p_Centers_thread[nCurThread][i][j];
					p_nLocalMembers[i] += p_nMembers_thread[nCurThread][i];
				}
		}
		
		// recalculating class centers in master
		MPI_Reduce(p_nLocalMembers, p_nMembers, m_nClass, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		
		for (int i=0; i<m_nClass; i++)
			MPI_Reduce(p_LocalCenters[i], p_Centers[i], m_nBand, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
		if (m_nID == 0)
		{
			for (int i=0; i<m_nClass; i++)
				for (int j=0; j<m_nBand; j++)
				{
					if (p_nMembers[i] != 0)
						p_Centers[i][j] = p_Centers[i][j]/p_nMembers[i];
					else
						p_Centers[i][j] = -1;

				}
		}

		// transferring class centers
		for (int i=0; i<m_nClass; i++)
			MPI_Bcast(p_Centers[i], m_nBand, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		MPI_Reduce(&nLocalChanged, &nTotalChanged, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

		gettimeofday(&t_end, 0);
		if (m_nID == 0)
		{
			fprintf(stderr, "%d th iteration: Classified in %.2f seconds and %.4lf changed\n",
				iter, float(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_usec-t_start.tv_usec)/1000000.0, double(nTotalChanged)/double(m_nEffectivePixel));
			fflush(stderr);
			fprintf(pf, "%d th iteration: Classified in %.2f seconds and %.4lf changed\n",
				iter, float(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_usec-t_start.tv_usec)/1000000.0, double(nTotalChanged)/double(m_nEffectivePixel));
		}
		
		// writing stats				
		if (m_nID == 0)
		{
			for (int i=0; i<m_nClass; i++)
				for (int j=0; j<m_nBand; j++)
					fprintf(pf, "class %d band %d center %.2lf member %d \n", i, j, p_Centers[i][j], p_nMembers[i]);
		}
				

		// terminating condition
		if (double(nTotalChanged)/double(m_nEffectivePixel)<m_nChangedRatio)
			m_bIteration = 0;

		MPI_Bcast(&m_bIteration, 1, MPI_INT, 0, MPI_COMM_WORLD);

		if (m_bIteration == 0)
			break;
	}


	// 5. writing results
	if (m_nID == 0)
	{
		unsigned char *p_Class;
		p_Class = new unsigned char[m_nPixelsPerNode];
		p_file = fopen(outPath, "wb");

		gettimeofday(&t_start, 0);
		for (int i=1; i<=m_nProcessor-2; i++)
		{
			MPI_Recv(p_Class, m_nPixelsPerNode, MPI_UNSIGNED_CHAR, i, 200, MPI_COMM_WORLD, &status);
			fwrite(p_Class, sizeof(unsigned char), m_nPixelsPerNode, p_file);
		}
		int nRemaining = m_nRow*m_nColumn-(m_nProcessor-2)*m_nPixelsPerNode;
		MPI_Recv(p_Class, nRemaining, MPI_UNSIGNED_CHAR, m_nProcessor-1, 200, MPI_COMM_WORLD, &status);
		fwrite(p_Class, sizeof(unsigned char), nRemaining, p_file);

		fclose(p_file);

		gettimeofday(&t_end, 0);
		fprintf(stderr, "Result is transferred and written in %.2f seconds\n", float(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_usec-t_start.tv_usec)/1000000.0);
		fflush(stderr);
		fprintf(pf, "Result is transferred and written in %.2f seconds\n", float(t_end.tv_sec-t_start.tv_sec)+(t_end.tv_usec-t_start.tv_usec)/1000000.0);
	}
	else
	{
		MPI_Send(p_Class, m_nPixelsInNode, MPI_UNSIGNED_CHAR, 0, 200, MPI_COMM_WORLD);
	}

		
	// 6. releasing arrays
	if (m_nID != 0)
	{
		for(int i=0; i<m_nBand; i++)
			delete p_Image[i];
		delete p_Image;

		for (int nCurThread = 0; nCurThread<m_nThread; nCurThread++)
		{
			for(int i=0; i<m_nClass; i++)
				delete p_Centers_thread[nCurThread][i];

			delete p_nMembers_thread[nCurThread];
		}
		delete p_Centers_thread;
		delete p_nMembers_thread;
	}

	for(int i=0; i<m_nClass; i++)
		delete p_Centers[i], p_LocalCenters[i];
	delete p_Centers, p_LocalCenters;

	delete p_nMembers, p_nLocalMembers;

	// 7. finalizing
	if (m_nID == 0)
	{
		gettimeofday(&t_end0, 0);
		fprintf(stderr, "The whole process is done in %.2f seconds\n", float(t_end0.tv_sec-t_start0.tv_sec)+(t_end0.tv_usec-t_start0.tv_usec)/1000000.0);
		fflush(stderr);
		fprintf(pf, "The whole process is done in %.2f seconds\n", float(t_end0.tv_sec-t_start0.tv_sec)+(t_end0.tv_usec-t_start0.tv_usec)/1000000.0);
	}
	
	if (m_nID == 0)
		fclose(pf);

	MPI_Finalize();

	return 0;
}

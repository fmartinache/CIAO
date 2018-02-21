/**
 * @file    ImageCreate.c
 * @brief   Create image
 * 
 * Creates images and streams (shared memory)
 *  
 * @author  O. Guyon
 * @date    10 Jul 2017
 *
 * 
 * @bug No known bugs.
 * 
 */



#define _GNU_SOURCE

#include <stdint.h>
#include <unistd.h>
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <sys/mman.h>
#include <signal.h> 

#include <semaphore.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <fcntl.h> // for open
#include <unistd.h> // for close
#include <errno.h>

#include "ImageStruct.h"




#ifdef __MACH__
#include <mach/mach_time.h>
#define CLOCK_REALTIME 0
#define CLOCK_MONOTONIC 0
static int clock_gettime(int clk_id, struct mach_timespec *t){
    mach_timebase_info_data_t timebase;
    mach_timebase_info(&timebase);
    uint64_t time;
    time = mach_absolute_time();
    double nseconds = ((double)time * (double)timebase.numer)/((double)timebase.denom);
    double seconds = ((double)time * (double)timebase.numer)/((double)timebase.denom * 1e9);
    t->tv_sec = seconds;
    t->tv_nsec = nseconds;
    return 0;
}
#else
#include <time.h>
#endif
















int ImageCreate_printERROR(const char *file, const char *func, int line, char *errmessage)
{
    char buff[256];

    fprintf(stderr,"%c[%d;%dm ERROR [ FILE: %s   FUNCTION: %s   LINE: %d ]  %c[%d;m\n", (char) 27, 1, 31, file, func, line, (char) 27, 0);
    if( errno != 0)
    {
        if( strerror_r( errno, buff, 256 ) == 0 ) {
            fprintf(stderr,"C Error: %s\n", buff );
        }
        else
            fprintf(stderr,"Unknown C Error\n");
    }
    else
        fprintf(stderr,"No C error (errno = 0)\n");

    fprintf(stderr,"%c[%d;%dm %s  %c[%d;m\n", (char) 27, 1, 31, errmessage, (char) 27, 0);

    return(0);
}









int ImageCreateSem(IMAGE *image, long NBsem)
{
    char sname[200];
    long s, s1;
    int r;
    char command[200];
    char fname[200];
    int semfile[100];

	
	// Remove pre-existing semaphores if any
    if(image->md[0].sem != NBsem)
    {
        // Close existing semaphores ...
        for(s=0; s < image->md[0].sem; s++)
            sem_close(image->semptr[s]);
        image->md[0].sem = 0;

		// ... and remove associated files
        for(s1=NBsem; s1<100; s1++)
        {
            sprintf(fname, "/dev/shm/sem.%s_sem%02ld", image->md[0].name, s1);
            remove(fname);
        }
        free(image->semptr);
        image->semptr = NULL;
    }

   
    if(image->md[0].sem == 0)
    {
        if(image->semptr!=NULL)
            free(image->semptr);

        image->md[0].sem = NBsem;
        printf("malloc semptr %d entries\n", image->md[0].sem);
        image->semptr = (sem_t**) malloc(sizeof(sem_t**)*image->md[0].sem);


        for(s=0; s<NBsem; s++)
        {
            sprintf(sname, "%s_sem%02ld", image->md[0].name, s);
            if ((image->semptr[s] = sem_open(sname, 0, 0644, 0))== SEM_FAILED) {
                if ((image->semptr[s] = sem_open(sname, O_CREAT, 0644, 1)) == SEM_FAILED) {
                    perror("semaphore initilization");
                }
                else
                    sem_init(image->semptr[s], 1, 0);
            }
        }
    }
    
    return(0);
}





int ImageCreate(IMAGE *image, const char *name, long naxis, uint32_t *size, uint8_t atype, int shared, int NBkw)
{
    long i,ii;
    time_t lt;
    long nelement;
    struct timespec timenow;
    char sname[200];

    size_t sharedsize = 0; // shared memory size in bytes
    int SM_fd; // shared memory file descriptor
    char SM_fname[200];
    int result;
    IMAGE_METADATA *map;
    char *mapv; // pointed cast in bytes

    int kw;
    char comment[80];
    char kname[16];
    

    nelement = 1;
    for(i=0; i<naxis; i++)
        nelement*=size[i];

    // compute total size to be allocated
    if(shared==1)
    {
        // create semlog

        sprintf(sname, "%s_semlog", name);
        remove(sname);
        image->semlog = NULL;

        if ((image->semlog = sem_open(sname, O_CREAT, 0644, 1)) == SEM_FAILED)
            perror("semaphore creation / initilization");
        else
            sem_init(image->semlog, 1, 0);



        sharedsize = sizeof(IMAGE_METADATA);

        if(atype == _DATATYPE_UINT8)
            sharedsize += nelement*SIZEOF_DATATYPE_UINT8;
        if(atype == _DATATYPE_INT8)
            sharedsize += nelement*SIZEOF_DATATYPE_INT8;

        if(atype == _DATATYPE_UINT16)
            sharedsize += nelement*SIZEOF_DATATYPE_UINT16;
        if(atype == _DATATYPE_INT16)
            sharedsize += nelement*SIZEOF_DATATYPE_INT16;

        if(atype == _DATATYPE_INT32)
            sharedsize += nelement*SIZEOF_DATATYPE_INT32;
        if(atype == _DATATYPE_UINT32)
            sharedsize += nelement*SIZEOF_DATATYPE_UINT32;


        if(atype == _DATATYPE_INT64)
            sharedsize += nelement*SIZEOF_DATATYPE_INT64;

        if(atype == _DATATYPE_UINT64)
            sharedsize += nelement*SIZEOF_DATATYPE_UINT64;


        if(atype == _DATATYPE_FLOAT)
            sharedsize += nelement*SIZEOF_DATATYPE_FLOAT;

        if(atype == _DATATYPE_DOUBLE)
            sharedsize += nelement*SIZEOF_DATATYPE_DOUBLE;

        if(atype == _DATATYPE_COMPLEX_FLOAT)
            sharedsize += nelement*SIZEOF_DATATYPE_COMPLEX_FLOAT;

        if(atype == _DATATYPE_COMPLEX_DOUBLE)
            sharedsize += nelement*SIZEOF_DATATYPE_COMPLEX_DOUBLE;


        sharedsize += NBkw*sizeof(IMAGE_KEYWORD);


        sprintf(SM_fname, "%s/%s.im.shm", SHAREDMEMDIR, name);
        SM_fd = open(SM_fname, O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);

        if (SM_fd == -1) {
            perror("Error opening file for writing");
            exit(0);
        }

        image->shmfd = SM_fd;
        image->memsize = sharedsize;

        result = lseek(SM_fd, sharedsize-1, SEEK_SET);
        if (result == -1) {
            close(SM_fd);
            ImageCreate_printERROR(__FILE__,__func__,__LINE__,"Error calling lseek() to 'stretch' the file");
            exit(0);
        }

        result = write(SM_fd, "", 1);
        if (result != 1) {
            close(SM_fd);
            perror("Error writing last byte of the file");
            exit(0);
        }

        map = (IMAGE_METADATA*) mmap(0, sharedsize, PROT_READ | PROT_WRITE, MAP_SHARED, SM_fd, 0);
        if (map == MAP_FAILED) {
            close(SM_fd);
            perror("Error mmapping the file");
            exit(0);
        }

        image->md = (IMAGE_METADATA*) map;
        image->md[0].sem = 0;
        image->md[0].shared = 1;
    }
    else
    {
        image->md = (IMAGE_METADATA*) malloc(sizeof(IMAGE_METADATA));
        image->md[0].shared = 0;
        if(NBkw>0)
            image->kw = (IMAGE_KEYWORD*) malloc(sizeof(IMAGE_KEYWORD)*NBkw);
        else
            image->kw = NULL;
    }

    image->md[0].atype = atype;
    image->md[0].naxis = naxis;
    strcpy(image->name, name); // local name
    strcpy(image->md[0].name, name);
    for(i=0; i<naxis; i++)
        image->md[0].size[i] = size[i];
    image->md[0].NBkw = NBkw;


    if(atype == _DATATYPE_UINT8)
    {
        if(shared==1)
        {
			mapv = (char*) map;
            mapv += sizeof(IMAGE_METADATA);
            image->array.UI8 = (uint8_t*) (mapv);
            memset(image->array.UI8, '\0', nelement*sizeof(uint8_t));
            mapv += sizeof(uint8_t)*nelement;
            image->kw = (IMAGE_KEYWORD*) (mapv);
		}
        else
            image->array.UI8 = (uint8_t*) calloc ((size_t) nelement, sizeof(uint8_t));


        if(image->array.UI8 == NULL)
        {
            ImageCreate_printERROR(__FILE__,__func__,__LINE__,"memory allocation failed");
            fprintf(stderr,"%c[%d;%dm", (char) 27, 1, 31);
            fprintf(stderr,"Image name = %s\n",name);
            fprintf(stderr,"Image size = ");
            fprintf(stderr,"%ld", (long) size[0]);
            for(i=1; i<naxis; i++)
                fprintf(stderr,"x%ld", (long) size[i]);
            fprintf(stderr,"\n");
            fprintf(stderr,"Requested memory size = %ld elements = %f Mb\n", (long) nelement,1.0/1024/1024*nelement*sizeof(uint8_t));
            fprintf(stderr," %c[%d;m",(char) 27, 0);
            exit(0);
        }
    }

    if(atype == _DATATYPE_INT8)
    {
        if(shared==1)
        {
			mapv = (char*) map;
            mapv += sizeof(IMAGE_METADATA);
            image->array.SI8 = (int8_t*) (mapv);
            memset(image->array.SI8, '\0', nelement*sizeof(int8_t));
            mapv += sizeof(int8_t)*nelement;
            image->kw = (IMAGE_KEYWORD*) (mapv);
		}
        else
            image->array.SI8 = (int8_t*) calloc ((size_t) nelement, sizeof(int8_t));


        if(image->array.SI8 == NULL)
        {
            ImageCreate_printERROR(__FILE__,__func__,__LINE__,"memory allocation failed");
            fprintf(stderr,"%c[%d;%dm", (char) 27, 1, 31);
            fprintf(stderr,"Image name = %s\n",name);
            fprintf(stderr,"Image size = ");
            fprintf(stderr,"%ld", (long) size[0]);
            for(i=1; i<naxis; i++)
                fprintf(stderr,"x%ld", (long) size[i]);
            fprintf(stderr,"\n");
            fprintf(stderr,"Requested memory size = %ld elements = %f Mb\n", (long) nelement,1.0/1024/1024*nelement*sizeof(int8_t));
            fprintf(stderr," %c[%d;m",(char) 27, 0);
            exit(0);
        }
    }



    if(atype == _DATATYPE_UINT16)
    {
        if(shared==1)
        {
            mapv = (char*) map;
            mapv += sizeof(IMAGE_METADATA);
            image->array.UI16 = (uint16_t*) (mapv);
            memset(image->array.UI16, '\0', nelement*sizeof(uint16_t));
            mapv += sizeof(uint16_t)*nelement;
            image->kw = (IMAGE_KEYWORD*) (mapv);
        }
        else
            image->array.UI16 = (uint16_t*) calloc ((size_t) nelement, sizeof(uint16_t));

        if(image->array.UI16 == NULL)
        {
            ImageCreate_printERROR(__FILE__,__func__,__LINE__,"memory allocation failed");
            fprintf(stderr,"%c[%d;%dm", (char) 27, 1, 31);
            fprintf(stderr,"Image name = %s\n",name);
            fprintf(stderr,"Image size = ");
            fprintf(stderr,"%ld", (long) size[0]);
            for(i=1; i<naxis; i++)
                fprintf(stderr,"x%ld", (long) size[i]);
            fprintf(stderr,"\n");
            fprintf(stderr,"Requested memory size = %ld elements = %f Mb\n", (long) nelement, 1.0/1024/1024*nelement*sizeof(uint16_t));
            fprintf(stderr," %c[%d;m",(char) 27, 0);
            exit(0);
        }
    }

    if(atype == _DATATYPE_INT16)
    {
        if(shared==1)
        {
            mapv = (char*) map;
            mapv += sizeof(IMAGE_METADATA);
            image->array.SI16 = (int16_t*) (mapv);
            memset(image->array.SI16, '\0', nelement*sizeof(int16_t));
            mapv += sizeof(int16_t)*nelement;
            image->kw = (IMAGE_KEYWORD*) (mapv);
        }
        else
            image->array.SI16 = (int16_t*) calloc ((size_t) nelement, sizeof(int16_t));

        if(image->array.SI16 == NULL)
        {
            ImageCreate_printERROR(__FILE__,__func__,__LINE__,"memory allocation failed");
            fprintf(stderr,"%c[%d;%dm", (char) 27, 1, 31);
            fprintf(stderr,"Image name = %s\n",name);
            fprintf(stderr,"Image size = ");
            fprintf(stderr,"%ld", (long) size[0]);
            for(i=1; i<naxis; i++)
                fprintf(stderr,"x%ld", (long) size[i]);
            fprintf(stderr,"\n");
            fprintf(stderr,"Requested memory size = %ld elements = %f Mb\n", (long) nelement, 1.0/1024/1024*nelement*sizeof(int16_t));
            fprintf(stderr," %c[%d;m",(char) 27, 0);
            exit(0);
        }
    }


    if(atype == _DATATYPE_UINT32)
    {
        if(shared==1)
        {
			mapv = (char*) map;
            mapv += sizeof(IMAGE_METADATA);
            image->array.UI32 = (uint32_t*) (mapv);
            memset(image->array.UI32, '\0', nelement*sizeof(uint32_t));
            mapv += sizeof(uint32_t)*nelement;
            image->kw = (IMAGE_KEYWORD*) (mapv);
        }
        else
            image->array.UI32 = (uint32_t*) calloc ((size_t) nelement, sizeof(uint32_t));

        if(image->array.UI32 == NULL)
        {
            ImageCreate_printERROR(__FILE__,__func__,__LINE__,"memory allocation failed");
            fprintf(stderr,"%c[%d;%dm", (char) 27, 1, 31);
            fprintf(stderr,"Image name = %s\n",name);
            fprintf(stderr,"Image size = ");
            fprintf(stderr,"%ld", (long) size[0]);
            for(i=1; i<naxis; i++)
                fprintf(stderr,"x%ld", (long) size[i]);
            fprintf(stderr,"\n");
            fprintf(stderr,"Requested memory size = %ld elements = %f Mb\n", (long) nelement,1.0/1024/1024*nelement*sizeof(uint32_t));
            fprintf(stderr," %c[%d;m",(char) 27, 0);
            exit(0);
        }
    }
    


    if(atype == _DATATYPE_INT32)
    {
        if(shared==1)
        {
			mapv = (char*) map;
            mapv += sizeof(IMAGE_METADATA);
            image->array.SI32 = (int32_t*) (mapv);
            memset(image->array.SI32, '\0', nelement*sizeof(int32_t));
            mapv += sizeof(int32_t)*nelement;
            image->kw = (IMAGE_KEYWORD*) (mapv);
        }
        else
            image->array.SI32 = (int32_t*) calloc ((size_t) nelement, sizeof(int32_t));

        if(image->array.SI32 == NULL)
        {
            ImageCreate_printERROR(__FILE__,__func__,__LINE__,"memory allocation failed");
            fprintf(stderr,"%c[%d;%dm", (char) 27, 1, 31);
            fprintf(stderr,"Image name = %s\n",name);
            fprintf(stderr,"Image size = ");
            fprintf(stderr,"%ld", (long) size[0]);
            for(i=1; i<naxis; i++)
                fprintf(stderr,"x%ld", (long) size[i]);
            fprintf(stderr,"\n");
            fprintf(stderr,"Requested memory size = %ld elements = %f Mb\n", (long) nelement,1.0/1024/1024*nelement*sizeof(int32_t));
            fprintf(stderr," %c[%d;m",(char) 27, 0);
            exit(0);
        }
    }
    
    
    
    if(atype == _DATATYPE_UINT64)
    {
        if(shared==1)
        {
			mapv = (char*) map;
            mapv += sizeof(IMAGE_METADATA);
            image->array.UI64 = (uint64_t*) (mapv);
            memset(image->array.UI64, '\0', nelement*sizeof(uint64_t));
            mapv += sizeof(uint64_t)*nelement;
            image->kw = (IMAGE_KEYWORD*) (mapv);
        }
        else
            image->array.UI64 = (uint64_t*) calloc ((size_t) nelement, sizeof(uint64_t));

        if(image->array.SI64 == NULL)
        {
            ImageCreate_printERROR(__FILE__,__func__,__LINE__,"memory allocation failed");
            fprintf(stderr,"%c[%d;%dm", (char) 27, 1, 31);
            fprintf(stderr,"Image name = %s\n",name);
            fprintf(stderr,"Image size = ");
            fprintf(stderr,"%ld", (long) size[0]);
            for(i=1; i<naxis; i++)
                fprintf(stderr,"x%ld", (long) size[i]);
            fprintf(stderr,"\n");
            fprintf(stderr,"Requested memory size = %ld elements = %f Mb\n", (long) nelement,1.0/1024/1024*nelement*sizeof(uint64_t));
            fprintf(stderr," %c[%d;m",(char) 27, 0);
            exit(0);
        }
    }
    
    if(atype == _DATATYPE_INT64)
    {
        if(shared==1)
        {
			mapv = (char*) map;
            mapv += sizeof(IMAGE_METADATA);
            image->array.SI64 = (int64_t*) (mapv);
            memset(image->array.SI64, '\0', nelement*sizeof(int64_t));
            mapv += sizeof(int64_t)*nelement;
            image->kw = (IMAGE_KEYWORD*) (mapv);
        }
        else
            image->array.SI64 = (int64_t*) calloc ((size_t) nelement, sizeof(int64_t));

        if(image->array.SI64 == NULL)
        {
            ImageCreate_printERROR(__FILE__,__func__,__LINE__,"memory allocation failed");
            fprintf(stderr,"%c[%d;%dm", (char) 27, 1, 31);
            fprintf(stderr,"Image name = %s\n",name);
            fprintf(stderr,"Image size = ");
            fprintf(stderr,"%ld", (long) size[0]);
            for(i=1; i<naxis; i++)
                fprintf(stderr,"x%ld", (long) size[i]);
            fprintf(stderr,"\n");
            fprintf(stderr,"Requested memory size = %ld elements = %f Mb\n", (long) nelement,1.0/1024/1024*nelement*sizeof(int64_t));
            fprintf(stderr," %c[%d;m",(char) 27, 0);
            exit(0);
        }
    }


    if(atype == _DATATYPE_FLOAT)	{
        if(shared==1)
        {
            mapv = (char*) map;
            mapv += sizeof(IMAGE_METADATA);
            image->array.F = (float*) (mapv);
            memset(image->array.F, '\0', nelement*sizeof(float));
            mapv += sizeof(float)*nelement;
            image->kw = (IMAGE_KEYWORD*) (mapv);
        }
        else
            image->array.F = (float*) calloc ((size_t) nelement, sizeof(float));

        if(image->array.F == NULL)
        {
            ImageCreate_printERROR(__FILE__,__func__,__LINE__,"memory allocation failed");
            fprintf(stderr,"%c[%d;%dm", (char) 27, 1, 31);
            fprintf(stderr,"Image name = %s\n",name);
            fprintf(stderr,"Image size = ");
            fprintf(stderr,"%ld", (long) size[0]);
            for(i=1; i<naxis; i++)
                fprintf(stderr,"x%ld", (long) size[i]);
            fprintf(stderr,"\n");
            fprintf(stderr,"Requested memory size = %ld elements = %f Mb\n", (long) nelement,1.0/1024/1024*nelement*sizeof(float));
            fprintf(stderr," %c[%d;m",(char) 27, 0);
            exit(0);
        }
    }
    
    if(atype == _DATATYPE_DOUBLE)
    {
        if(shared==1)
        {
			mapv = (char*) map;
			mapv += sizeof(IMAGE_METADATA);
			image->array.D = (double*) (mapv);
            memset(image->array.D, '\0', nelement*sizeof(double));
            mapv += sizeof(double)*nelement;
            image->kw = (IMAGE_KEYWORD*) (mapv);
        }
        else
            image->array.D = (double*) calloc ((size_t) nelement, sizeof(double));
  
        if(image->array.D == NULL)
        {
            ImageCreate_printERROR(__FILE__,__func__,__LINE__,"memory allocation failed");
            fprintf(stderr,"%c[%d;%dm", (char) 27, 1, 31);
            fprintf(stderr,"Image name = %s\n",name);
            fprintf(stderr,"Image size = ");
            fprintf(stderr,"%ld", (long) size[0]);
            for(i=1; i<naxis; i++)
                fprintf(stderr,"x%ld", (long) size[i]);
            fprintf(stderr,"\n");
            fprintf(stderr,"Requested memory size = %ld elements = %f Mb\n", (long) nelement,1.0/1024/1024*nelement*sizeof(double));
            fprintf(stderr," %c[%d;m",(char) 27, 0);
            exit(0);
        }
    }
    
    if(atype == _DATATYPE_COMPLEX_FLOAT)
    {
        if(shared==1)
        {
			mapv = (char*) map;
			mapv += sizeof(IMAGE_METADATA);
			image->array.CF = (complex_float*) (mapv);			
            memset(image->array.CF, '\0', nelement*sizeof(complex_float));
            mapv += sizeof(complex_float)*nelement;
            image->kw = (IMAGE_KEYWORD*) (mapv);
        }
        else
            image->array.CF = (complex_float*) calloc ((size_t) nelement, sizeof(complex_float));

        if(image->array.CF == NULL)
        {
            ImageCreate_printERROR(__FILE__,__func__,__LINE__,"memory allocation failed");
            fprintf(stderr,"%c[%d;%dm", (char) 27, 1, 31);
            fprintf(stderr,"Image name = %s\n",name);
            fprintf(stderr,"Image size = ");
            fprintf(stderr,"%ld", (long) size[0]);
            for(i=1; i<naxis; i++)
                fprintf(stderr,"x%ld", (long) size[i]);
            fprintf(stderr,"\n");
            fprintf(stderr,"Requested memory size = %ld elements = %f Mb\n", (long) nelement,1.0/1024/1024*nelement*sizeof(complex_float));
            fprintf(stderr," %c[%d;m",(char) 27, 0);
            exit(0);
        }
    }
    
    if(atype == _DATATYPE_COMPLEX_DOUBLE)
    {
        if(shared==1)
        {
			mapv = (char*) map;
			mapv += sizeof(IMAGE_METADATA);
			image->array.CD = (complex_double*) (mapv);			
            memset(image->array.CD, '\0', nelement*sizeof(complex_double));
            mapv += sizeof(complex_double)*nelement;
            image->kw = (IMAGE_KEYWORD*) (mapv);
        }
        else
            image->array.CD = (complex_double*) calloc ((size_t) nelement,sizeof(complex_double));

        if(image->array.CD == NULL)
        {
            ImageCreate_printERROR(__FILE__,__func__,__LINE__,"memory allocation failed");
            fprintf(stderr,"%c[%d;%dm", (char) 27, 1, 31);
            fprintf(stderr,"Image name = %s\n",name);
            fprintf(stderr,"Image size = ");
            fprintf(stderr,"%ld", (long) size[0]);
            for(i=1; i<naxis; i++)
                fprintf(stderr,"x%ld", (long) size[i]);
            fprintf(stderr,"\n");
            fprintf(stderr,"Requested memory size = %ld elements = %f Mb\n", (long) nelement,1.0/1024/1024*nelement*sizeof(complex_double));
            fprintf(stderr," %c[%d;m",(char) 27, 0);
            exit(0);
        }
    }



    clock_gettime(CLOCK_REALTIME, &timenow);
    image->md[0].last_access = 1.0*timenow.tv_sec + 0.000000001*timenow.tv_nsec;
    image->md[0].creation_time = image->md[0].last_access;
    image->md[0].write = 0;
    image->md[0].cnt0 = 0;
    image->md[0].cnt1 = 0;
    image->md[0].nelement = nelement;

    if(shared==1)
        ImageCreateSem(image, 10); // by default, create 10 semaphores


    // initialize keywords
    for(kw=0; kw<image->md[0].NBkw; kw++)
        image->kw[kw].type = 'N';


    return(0);
}


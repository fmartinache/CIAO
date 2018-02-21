/**
 * @file    ImageCreate.h
 * @brief   Function prototypes for ImageCreate
 * 
 *  
 * @author  O. Guyon
 * @date    9 Jul 2017
 *
 * 
 * @bug No known bugs.
 * 
 */



#ifndef _IMAGECREATE_H
#define _IMAGECREATE_H




int ImageCreateSem(IMAGE *image, long NBsem);

int ImageCreate(IMAGE *image, const char *name, long naxis, uint32_t *size, uint8_t atype, int shared, int NBkw);


#endif

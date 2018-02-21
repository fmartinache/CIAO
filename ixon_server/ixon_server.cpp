/* =========================================================================
   Source code for the camera server running the andor iXon camera for CIAO.
   The server awaits from simple commands like "stream", "quit", "tint xxx"
   on a named pipe (aka a fifo, for "first in, first out") and executes them
   assuming that the command is valid.
   
   The goal of this version is to implement the revised shared memory data
   structure developed during the summer 2017.

   Frantz Martinache.
   ========================================================================= */

#include <stdio.h>

#ifdef __GNUC__
#  if(__GNUC__ > 3 || __GNUC__ ==3)
#	define _GNUC3_
#  endif
#endif

#ifdef _GNUC3_
#  include <iostream>
#  include <fstream>
   using namespace std;
#else
#  include <iostream.h>
#  include <fstream.h>
#endif

#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>

#include "atmcdLXd.h"
#include <fitsio.h>

//#include "SHM_DS.h"

#include "ImageStruct.h"
#include "ImageCreate.h"

#include <pthread.h>
#include <time.h>

// ==========================================================================
//                              CONSTANTS 
// ==========================================================================
#define LINESZ 256                         // string length
#define IXON_SM_FNAME "/tmp/ixon.shm"      // shared memory file for iXon
#define msleep(x) usleep((int)(1000*(x))); // sleep in milli-seconds

// ==========================================================================
//                  structure storing acquisition information
// ==========================================================================
typedef struct {
  float exp_time; // exposure time in seconds
  float acc_time; // accumulate cycle time in seconds (cf SDK)
  float kin_time; // kinetic cycle time in seconds (cf SDK)
  // ---------
  int XW;         // image (window) width
  int YW;         // image (window) height
  // ---------
  int nleft;      // number of images left in acquisition
  int nextframe;  // index of the next frame to write
  // ---------
  bool bstreamon; // camera streaming? (continuous acq.)
  bool bacqon;    // camera acquiring? 
  bool babort;    // abort command was issued?
  bool bcamOK;    // the happy camera flag!
} cam_config;

// ==========================================================================
//                          GLOBAL VARIABLES
// ==========================================================================
cam_config *cconf;                       // pointer to acq. info structure

char fits_name[LINESZ] = "test.fits";    // temporary feature for tests
char conf_name[LINESZ] = "config.xml";   // OptAcqr config file

char myfifoin[LINESZ] = "/home/ciaodev/bin/ixon_fifo_in";  // out-going pipe
char myfifout[LINESZ] = "/home/ciaodev/bin/ixon_fifo_out"; // out-going pipe

int first_im = 0;
int last_im = 0;

// size_t shm_size = 0;                     // shared memory size in bytes
// int fd_shm = -1;                         // file descriptor for shared memory
// SHM_DATA *shm_data;                      // pointer to shared memory for iXon

IMAGE *imarray;            // pointer to image

// ==========================================================================
//                              PROTOTYPES
// ==========================================================================
int CameraSelect(int iNumArgs, char* szArgList[]);
void cam_config_init(cam_config *camconf);
void printBits(unsigned int num);

// ==========================================================================
//                         FUNCTIONS 
// ==========================================================================
void printBits(unsigned int num) {
  for(unsigned int bit=0;bit<(sizeof(unsigned int) * 8); bit++) {
    printf("%i ", num & 0x01);
    num = num >> 1;
  }
  printf("\n");
}

// ---------------------------------------------------------
//        initialize the camconf structure
// ---------------------------------------------------------
void cam_config_init(cam_config* camconf) {
  camconf->exp_time = 0.00001; // default: shortest exp time
  camconf->kin_time = 0.0;
  camconf->bstreamon = false;
  camconf->bacqon = false;
  camconf->babort = false;
  camconf->bcamOK = true;
}

// ==========================================================================
//                         ACQUISITION THREAD
// ==========================================================================
void* acquire(void *params) { // camera acquisition thread
  cam_config* camconf = (cam_config*) params;
  int status;
  int nel = camconf->XW * camconf->YW;
  int nfr = 10, ifr = 0;     // variables used to estimate frame rate
  float t0 = 0.0, t1 = 0.0;  // time variables for frame rate
  int idisp = 0, ndisp = 20; // control the text output refresh rate
  float ifrate;              // frame rate (average over nfr frames)
  struct timespec now;       // clock readout
  struct tm* ptm;
  float *timing = (float*) malloc(nfr * sizeof(float));

  for (ifr = 0; ifr < nfr; ifr++) timing[ifr] = 0.1; // init timing array
  ifr = 1;

  // acquisition loop
  printf("\n");

  if (camconf->bstreamon) {
    StartAcquisition();
  }
  
  while (camconf->nleft > 0) {

    clock_gettime(CLOCK_REALTIME, &now);
    ptm = gmtime(&(now.tv_sec));
    t1 = (float)(now.tv_nsec)*1e-9 + (float)(ptm->tm_sec);

    // estimate the frame rate
    timing[ifr] = t1-t0;
    t0 = t1;
    ifr++;
    if (ifr == nfr) ifr = 0;

    ifrate = 0.0;
    for (int i = 0; i < nfr; i++)
      ifrate += timing[i];
    ifrate = (float)(nfr) / ifrate;

    //Loop until acquisition finished
    GetStatus(&status);
    //while(status==DRV_ACQUIRING) GetStatus(&status);

    // write to shared memory
    imarray[0].md[0].write = 1;                  // signal you are about to write
    GetMostRecentImage(imarray[0].array.SI32, nel); // direct write to SHM
    imarray[0].md[0].write = 0;                  // signal done writing data
    imarray[0].md[0].cnt0 ++;                    // increment counter!
    imarray[0].md[0].cnt1 ++;                    // increment counter!

    // check for abort
    if (camconf->babort) {
      camconf->nleft = 0;
      camconf->bstreamon = false;
      camconf->babort = false;
    }

    idisp++;
    if (idisp == ndisp) {
      printf("\r%02d Image [%9ld] frame rate = %.3f Hz, status code = %d", 
	     idisp, imarray[0].md[0].cnt0, ifrate, status);
      fflush(stdout);
      idisp = 0;
    }

    //GetNumberNewImages(&first_im, &last_im);
    //printf("  indices = %d %d", first_im, last_im);

    msleep(5); // sleep in milli-seconds
    if (!camconf->bstreamon) camconf->nleft--; // decrement if NOT streaming
  }
  printf("\n");

  camconf->bacqon    = false; // updating acquisition flags
  camconf->bstreamon = false; // updating acquisition flags
  return NULL;
}

// ==========================================================================
//                              MAIN PROGRAM
// ==========================================================================
int main(int argc, char* argv[]) {
  pthread_t tid_acqr;      // thread id for acquisition

  int mytest = -1;         // useful for ... tests
  long naxis = 2;          // number of axes
  uint8_t atype;           // data type code
  uint32_t *imsize;        // image size
  int shared;              // set to 1 if image is in shared memory
  int NBkw;                // number of keywords supported

  AndorCapabilities *caps; // TEST

  unsigned long error;
  bool quit;
  float fChoice;
  int gain = -1;
  int temp;
  int temp_min = 1, temp_max = 1;
  unsigned int state;

  int glow = -1, ghigh = -1;
  int nb_pre_gains = -1;
  int nb_adcs = -1, nb_amps = -1, nb_speeds = -1, 
    adc_channel = -1, amp_type = -1;
  float pre_gain = 0.0, hzs_speed = 0.0;

  char *cmd_in = (char*) malloc(LINESZ * sizeof(char));
  char *auxstr = (char*) malloc(LINESZ * sizeof(char));
  mytest *= 1;

  cconf = (cam_config*) malloc(sizeof(cam_config));
  cam_config_init(cconf);


  // ------------------------------------------
  //            Initialize CCD
  // ------------------------------------------
  if (CameraSelect (argc, argv) < 0) {
    cout << "*** CAMERA SELECTION ERROR" << endl;
    return -1;
  }
  
  error = Initialize((char*)"/usr/local/etc/andor/");
  if(error != DRV_SUCCESS){
    cout << "Initialisation error...exiting" << endl;
    return(1);
  }
  
  // ------------------------------------------
  //    startup configuration of the camera
  // ------------------------------------------
  sleep(2);                         // sleep to allow initialization to complete
  SetReadMode(4);                   // Set Read Mode to --Image--
  SetAcquisitionMode(5);            // Set Acq. mode to --Single scan--
  SetShutter(1,1,50,50);            // Initialize Shutter
  SetExposureTime(cconf->exp_time);           // Set initial exposure time
  GetDetector(&(cconf->XW), &(cconf->YW));    // Get Detector dimensions  
  SetImage(1, 1, 1, cconf->XW, 1, cconf->YW); // Setup Image dimensions

  SetKineticCycleTime(cconf->kin_time); // test

  error = GetAcquisitionTimings(&(cconf->exp_time), 
				&(cconf->acc_time), 
				&(cconf->kin_time));
  printf("\n >> Exp. time = %.5f s, acc. time = %.5f s, kin. time = %.5f s\n",
	 cconf->exp_time, cconf->acc_time, cconf->kin_time);

  caps = (AndorCapabilities *) malloc(sizeof(AndorCapabilities));
  caps->ulSize = sizeof(AndorCapabilities);
  error = GetCapabilities(caps);
  if(error != DRV_SUCCESS)
    cout << "GET Capabilities fail" << endl;

  cout << "Capabilities:" << endl << "============" << endl;
  printf("Acq Modes  : "); printBits(caps->ulAcqModes);
  printf("Read Modes : "); printBits(caps->ulReadModes);
  printf("Trig Modes : "); printBits(caps->ulTriggerModes);
  printf("Cam Type   : %u \n", caps->ulCameraType);
  printf("Pixl Mode  : "); printBits(caps->ulPixelMode);
  printf("Set Funcs  : "); printBits(caps->ulSetFunctions);
  printf("Get Funcs  : "); printBits(caps->ulGetFunctions);
  printf("Get Feats  : "); printBits(caps->ulFeatures);
  printf("PCI speed  : %u \n", caps->ulPCICard);
  printf("EMG capab  : "); printBits(caps->ulEMGainCapability);
  printf("FT Readmod : "); printBits(caps->ulFTReadModes);
  cout << "============" << endl;

  error = GetNumberPreAmpGains(&nb_pre_gains);
  printf("# of Preamp Gains : %d\n", nb_pre_gains);

  for (int i = 0; i < nb_pre_gains; i++) {
    error = GetPreAmpGain(i, &pre_gain);
    printf("%d: gain = %.2f\n", i, pre_gain);
  }
  error = SetPreAmpGain(0); // set pre-amp gain to 1.0

  // ====================================================================
  error = GetNumberAmp(&nb_amps);
  printf("%d output amplifier available \n", nb_amps);

  // ====================================================================
  error = GetNumberADChannels(&nb_adcs);
  printf("\n%d A/D channels are available \n", nb_adcs);

  int depth = -1;
  for (int i = 0; i < nb_adcs; i++) {
    error = GetBitDepth(i, &depth);
    printf("- channel #%d: depth = %d\n", i, depth);
  }
  printf("\n\n");
  // ====================================================================

  amp_type    = 0; // type of output amplification (1: conventional, 0: EMCCD)
  adc_channel = 1; //

  error = SetOutputAmplifier(amp_type);
  error = SetADChannel(adc_channel);

  printf("For AD channel = %d: ", adc_channel);
  error = GetNumberHSSpeeds(adc_channel, amp_type, &nb_speeds);
  printf("%d speeds are available in this mode\n", nb_speeds);

  for (int i = 0; i < nb_speeds; i++) {
    error = GetHSSpeed(adc_channel, amp_type, i, &hzs_speed);
    printf("- i=%d: speed = %.1f MHz\n", i, hzs_speed);
  }
  printf("\n");

  // ====================================================================

  //SetHSSpeed(amp_type, 0);
  //  SetOutputAmplifier(amp_type);

  //SetHighCapacity(1); // 
  error = GetEMCCDGain(&gain);

  if(error != DRV_SUCCESS)
    cout << "GET EMCCDGain fail" << endl;

  printf("EMCCD Gain = %d\n", gain);

  GetMCPGain(&gain);
  printf("MCP   Gain = %d\n", gain);

  GetMCPGainRange(&glow, &ghigh);
  printf("MCP Gain range = %d - %d\n", glow, ghigh);

  // ------------------------------------------
  //    special setup for the temperature
  // ------------------------------------------
  state = CoolerON();
  cout << "CoolerON status code: " << state << endl;

  state = GetTemperatureRange(&temp_min, &temp_max);
  cout << "Detector temperature range = " << temp_min << " - " << temp_max << endl;

  state = GetTemperature(&temp);    // get detector temperature
  cout << "Detector temperature = " << temp << " deg" << endl;

  temp = -5; // temperature setting (air cooling)
  state = SetTemperature(temp);

  state = GetTemperature(&temp);    // get detector temperature
  cout << "Detector temperature = " << temp << " deg" << endl;

  cout << "GetTemperature status code: " << state << endl;

  while(state != DRV_TEMPERATURE_STABILIZED) {
    sleep(1);
    state = GetTemperature(&temp);    // get detector temperature
    printf("\rTemperature = %+03d deg - status code: %6d", temp, state);
    fflush(stdout);
  }

  switch(state) {
  case DRV_TEMPERATURE_OFF: 
    cout << "Cooler OFF" << endl; break;
  case DRV_TEMPERATURE_STABILIZED: 
    cout << "Stabilized" << endl; break;
  case DRV_TEMPERATURE_NOT_REACHED: 
    cout << "Cooling" << endl; break;
  default:
    cout << "Cooler status unknown" << endl;
  }

  // ------------------------------------------------
  // create pipes for interaction with other programs
  // ------------------------------------------------
  if (mkfifo(myfifoin, 0777) != 0) printf("Could not create fifo!\n");
  if (mkfifo(myfifout, 0777) != 0) printf("Could not create fifo!\n");


  // ------------------------------------------
  //          setup shared memory 
  // ------------------------------------------
  naxis     = 2;                // 2D image
  imarray   = (IMAGE*) malloc(sizeof(IMAGE));
  imsize    = (uint32_t *) malloc(naxis * sizeof(uint32_t));
  imsize[0] = cconf->XW;       // image dimensions
  imsize[1] = cconf->YW;       // image dimensions
  atype     = _DATATYPE_INT32; // camera SDK writes "int"
  shared    = 1;               // image will be in shared memory
  NBkw      = 10;              // allocate space for 10 keywords
  ImageCreate(&imarray[0], "ixon", naxis, imsize, atype, shared, NBkw);


  // ------------------------------------------
  //          main server loop
  // ------------------------------------------
  int rfd = 0;
  int foo = -1;

  quit = false;
  printf("iXon server ready.");
  fflush(stdout);

  do {
    // ==============================================
    //    listen for commands on the input fifo
    // ==============================================
    rfd = open(myfifoin, O_RDONLY);
    sprintf(cmd_in, "%s", "          ");
    read(rfd, cmd_in, LINESZ);
    close(rfd);
    rfd = 0;

    // ===============================================
    //       processing the command from pipe
    // ===============================================

    // -----------------------------------
    //      abort current acquisition
    // -----------------------------------
    if (strncmp(cmd_in, "abort", 5) == 0) {
      if (cconf->bacqon) {
	cconf->babort = true;
	AbortAcquisition(); // camera command
      }
    }

    // -----------------------------------
    //         close the program
    // -----------------------------------
    if (strncmp(cmd_in, "quit", 4) == 0) {
      if (cconf->bacqon) {
	cconf->babort = true;
	AbortAcquisition(); // camera command
	printf("\nQuitting program\n");
	sleep(5);
      } 
      quit = true;
    }

    // -----------------------------------
    //       take a single image
    // -----------------------------------
    if (strncmp(cmd_in, "take ", 5) == 0) {
      foo = -1;
      sscanf(cmd_in, "%s %d", auxstr, &foo);
      cconf->nleft = (foo < 0) ? 1: foo; 
      cconf->bacqon = true;
      pthread_create(&tid_acqr, NULL, acquire, cconf);
    }

    // -----------------------------------
    //       continuous acquisition
    // -----------------------------------
    if (strncmp(cmd_in, "stream", 6) == 0) {
      if (!cconf->bstreamon) {
	cconf->bacqon = true;
	cconf->bstreamon = true;
	cconf->nleft = 1;
	pthread_create(&tid_acqr, NULL, acquire, cconf);
      }
    }

    // -----------------------------------
    //       update exposure time
    // -----------------------------------
    if (strncmp(cmd_in, "tint ", 5) == 0) {
      sscanf(cmd_in, "%s %f", auxstr, &fChoice);
      //printf("Exposure time: %f\n", fChoice);
      AbortAcquisition(); // camera command

      cconf->exp_time = fChoice;
      cconf->kin_time = 0.0;
      SetExposureTime(cconf->exp_time);
      SetKineticCycleTime(cconf->kin_time);

      error = GetAcquisitionTimings(&(cconf->exp_time), 
				    &(cconf->acc_time), 
				    &(cconf->kin_time));

      printf("Exp. time = %.5f s, acc. time = %.5f s, kin. time = %.5f s\n",
	     cconf->exp_time, cconf->acc_time, cconf->kin_time);

      sleep(1);
    }
  } while(!quit);
  
  printf("Shutting things down!\n");

  unlink(myfifoin);
  unlink(myfifout);

  fflush(stdout);

  free(imsize);
  free(imarray);

  ShutDown();                // Shut down the iXon
  free(cconf); cconf = NULL; // free cam_config structure
  return(EXIT_SUCCESS);
}

// --------------------------------------------------------------------------
//                              FUNCTION DEFINITIONS
// --------------------------------------------------------------------------
int CameraSelect (int iNumArgs, char* szArgList[]) {
  if (iNumArgs == 2) {
    
    at_32 lNumCameras;
    GetAvailableCameras(&lNumCameras);
    int iSelectedCamera = atoi(szArgList[1]);
    
    if (iSelectedCamera < lNumCameras && iSelectedCamera >= 0) {
      at_32 lCameraHandle;
      GetCameraHandle(iSelectedCamera, &lCameraHandle);
      SetCurrentCamera(lCameraHandle);
      return iSelectedCamera;
    }
    else
      return -1;
  }
  return(EXIT_SUCCESS);
}

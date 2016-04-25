/***************************************************
 ** ViSUS Visualization Project                    **
 ** Copyright (c) 2010 University of Utah          **
 ** Scientific Computing and Imaging Institute     **
 ** 72 S Central Campus Drive, Room 3750           **
 ** Salt Lake City, UT 84112                       **
 **                                                **
 ** For information about this project see:        **
 ** http://www.pascucci.org/visus/                 **
 **                                                **
 **      or contact: pascucci@sci.utah.edu         **
 **                                                **
 ****************************************************/

#include "PIDX_metadata.h"
#include <stdlib.h>
#include <tinyxml.h>

char* itoa(int val, int base){
  static char buf[32] = {0};
  int i = 30;
  for(; val && i ; --i, val /= base)
    buf[i] = "0123456789abcdef"[val % base];
  
  return &buf[i+1];
}

int PIDX_metadata_load(PIDX_metadata* metadata, const char* filename){
  *metadata = (PIDX_metadata)malloc(sizeof *(*metadata));
  memset(*metadata, 0, sizeof *(*metadata));
  strcpy((*metadata)->filename, filename);
  
  (*metadata)->doc = new TiXmlDocument((*metadata)->filename);
  bool loadOkay = (*metadata)->doc->LoadFile();
  if (loadOkay) { printf("XML loaded %s\n", (*metadata)->filename); return 0; }
  else          { printf("Could not load XML %s\n", (*metadata)->filename); return -1;}
  
  
}

int PIDX_metadata_create(PIDX_metadata* metadata, const char* filename){
  *metadata = (PIDX_metadata)malloc(sizeof *(*metadata));
  memset(*metadata, 0, sizeof *(*metadata));
  strcpy((*metadata)->filename, filename);
  
  (*metadata)->doc = new TiXmlDocument();
  TiXmlDeclaration * decl = new TiXmlDeclaration( "1.0", "", "" );
  (*metadata)->doc->LinkEndChild( decl );
  TiXmlElement * element = new TiXmlElement( "idx" );
  (*metadata)->doc->LinkEndChild( element );
  
  bool loadOkay = (*metadata)->doc->LoadFile();
  if (loadOkay) { printf("XML loaded %s\n", (*metadata)->filename); return 0; }
  else          { printf("Could not load XML %s\n", (*metadata)->filename); return -1;}
  
  return 0;
}

int PIDX_metadata_save(PIDX_metadata metadata){
  if(metadata != NULL){
    metadata->doc->SaveFile(metadata->filename);
    printf("file saved\n");
    return 0;
  }
  
  return -1;
};

int PIDX_metadata_add_timestep(PIDX_metadata metadata, int index, float value){
  if(!metadata)
    return -1;// use PIDX error
  
  TiXmlHandle docHandle( metadata->doc );
  
  TiXmlElement* timesteps = docHandle.FirstChild("idx").FirstChild("timesteps").ToElement();
  
  if(timesteps == NULL){
    TiXmlElement * root = metadata->doc->RootElement();
    
    timesteps = new TiXmlElement( "timesteps" );
    root->LinkEndChild( timesteps );
  }
  
  TiXmlElement * timestep = new TiXmlElement( "timestep" );
  timesteps->LinkEndChild(timestep);
  
  TiXmlText * log_time = new TiXmlText( itoa(index,10) );
  timestep->LinkEndChild(log_time);
  timestep->SetDoubleAttribute("time", value);
  
  return 0;
};

int PIDX_metadata_get_timestep(PIDX_metadata metadata, int index, float& value){
  if(!metadata)
    return -1;// use PIDX error
  
  TiXmlHandle docHandle( metadata->doc );
  
  TiXmlElement* child = docHandle.FirstChild("idx").FirstChild("timesteps").FirstChild("timestep").ToElement();
  
  for( ; child; child=child->NextSiblingElement() ){
    
    int time_log = atoi(child->GetText());
    
    if (time_log == index){
      value = atof(child->Attribute("time"));
      return 0;
    }
    
  }
  
  return -1;
}

int PIDX_metadata_destroy(PIDX_metadata metadata){
  if (metadata->doc != NULL) delete metadata->doc;
  if (metadata){ delete metadata; metadata = NULL; }
}



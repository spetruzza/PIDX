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

char* itoa(int val, int base)
{
  static char buf[32] = {0};
  if(val == 0){
    strcpy(buf,"0");
    return buf;
  }

  int i = 30;
  for(; val && i ; --i, val /= base)
    buf[i] = "0123456789abcdef"[val % base];
  
  return &buf[i+1];
}

PIDX_return_code PIDX_metadata_load(PIDX_metadata* metadata, const char* filename)
{
  *metadata = (PIDX_metadata)malloc(sizeof *(*metadata));
  memset(*metadata, 0, sizeof *(*metadata));
  strcpy((*metadata)->filename, filename);
  
  (*metadata)->doc = new TiXmlDocument((*metadata)->filename);
  bool loadOkay = (*metadata)->doc->LoadFile();
  if (loadOkay) { /*printf("XML loaded %s\n", (*metadata)->filename);*/ return PIDX_success; }
  else          { fprintf(stderr, "Could not load XML %s\n", (*metadata)->filename); return PIDX_err_metadata;}
  
}

PIDX_return_code PIDX_metadata_create(PIDX_metadata* metadata, const char* filename)
{
  *metadata = (PIDX_metadata)malloc(sizeof *(*metadata));
  memset(*metadata, 0, sizeof *(*metadata));
  strcpy((*metadata)->filename, filename);
  
  (*metadata)->doc = new TiXmlDocument();
  TiXmlDeclaration * decl = new TiXmlDeclaration( "1.0", "", "" );
  (*metadata)->doc->LinkEndChild( decl );
  TiXmlElement * element = new TiXmlElement( "idx" );
  (*metadata)->doc->LinkEndChild( element );
  
  bool loadOkay = (*metadata)->doc->LoadFile();
  if (loadOkay) { /*printf("XML loaded %s\n", (*metadata)->filename);*/ return PIDX_success; }
  else          { fprintf(stderr,"Could not load XML %s\n", (*metadata)->filename); return PIDX_err_metadata;}
}

PIDX_return_code PIDX_metadata_save(PIDX_metadata metadata)
{
  if(metadata != NULL){
    metadata->doc->SaveFile(metadata->filename);
    printf("file saved\n");
    return PIDX_success;
  }
  
  return PIDX_err_metadata;
};

PIDX_return_code PIDX_metadata_add_timestep(PIDX_metadata metadata, int index, double value)
{
  if(!metadata)
    return PIDX_err_metadata;// use PIDX error
  
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
  char phy_value[32];
  sprintf(phy_value,"%g",value);
  timestep->SetAttribute("time", phy_value);
    
  return PIDX_success;
};

PIDX_return_code PIDX_metadata_get_timestep(PIDX_metadata metadata, int index, double& value)
{
  if(!metadata)
    return PIDX_err_metadata;
  
  TiXmlHandle docHandle( metadata->doc );
  
  TiXmlElement* child = docHandle.FirstChild("idx").FirstChild("timesteps").FirstChild("timestep").ToElement();
  
  for( ; child; child=child->NextSiblingElement() )
  {  
    int time_log = atoi(child->GetText());
    
    if (time_log == index){
      value = strtod(child->Attribute("time"),NULL);
      return PIDX_success;
    }
    
  }
  
  return PIDX_err_metadata;
}

PIDX_return_code PIDX_metadata_add_simple_box(PIDX_metadata metadata, unsigned long long* log_size, float* phy_size)
{
  if(!metadata)
    return PIDX_err_metadata;// use PIDX error
  
  TiXmlHandle docHandle( metadata->doc );
  
  TiXmlElement* level = docHandle.FirstChild("idx").FirstChild("level").ToElement();
  
  if(level == NULL)
  {
    TiXmlElement * root = metadata->doc->RootElement();
    
    level = new TiXmlElement( "level" );
    root->LinkEndChild( level );
  }
  
  TiXmlElement * box = new TiXmlElement( "box" );
  level->LinkEndChild(box);

  TiXmlElement * lower = new TiXmlElement( "lower" );
  box->LinkEndChild(lower);
  TiXmlText * lower_value = new TiXmlText( "[0 0 0]" );
  lower->LinkEndChild(lower_value);

  TiXmlElement * upper = new TiXmlElement( "upper" );
  box->LinkEndChild(upper);
  char phy_values[256];
  sprintf(phy_values,"[%.5f %.5f %.5f]", phy_size[0], phy_size[1], phy_size[2]);
  TiXmlText * upper_value = new TiXmlText( phy_values );
  upper->LinkEndChild(upper_value);

  TiXmlElement * resolution = new TiXmlElement( "resolution" );
  box->LinkEndChild(resolution);
  char log_values[256];
  sprintf(log_values,"[%llu %llu %llu]", log_size[0], log_size[1], log_size[2]);
  TiXmlText * res_value = new TiXmlText( log_values );
  resolution->LinkEndChild(res_value);
    
  return PIDX_success;
}

PIDX_return_code PIDX_metadata_destroy(PIDX_metadata metadata)
{
  if (metadata->doc != NULL) delete metadata->doc;
  if (metadata){ delete metadata; metadata = NULL; }
  return PIDX_success;
}



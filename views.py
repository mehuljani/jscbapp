from django.shortcuts import render
from django.http import HttpResponse
from django.template import loader
from django.views.generic import ListView, CreateView, UpdateView
from django.urls import reverse_lazy
from django.http import HttpResponseRedirect
from django.http import StreamingHttpResponse
from django.conf import settings
from django.apps import apps
from django.core.mail import send_mail
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.views import View
from django.http import FileResponse, Http404
from datetime import datetime
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
from django.shortcuts import redirect

from .forms import *

import json
import os,csv
import pandas as pd
import sqlite3
import numpy as np
import json
import re
import subprocess as sub
import requests
import glob
import time
import base64

# Create your views here.

def file_upload(request):
    basedir = '/home/mehul/Mehul/djaenv/mysite/jscb/data/'
    if request.method == 'POST':
        form = FileUploadForm(request.POST, request.FILES)
        if form.is_valid():
            uploaded_file = form.cleaned_data['file']
            with open(basedir + uploaded_file.name, 'wb+') as destination:
                for chunk in uploaded_file.chunks():
                    destination.write(chunk)
            return redirect('file_list')
    else:
        form = FileUploadForm()
    return render(request, 'jscb/file_upload.html', {'form': form})	
    

@csrf_exempt
def results(request,selected_value=""):
    def imgProcess(imgpath):
        cmd = '''echo -n| base64 %s'''%(imgpath)
        out = str(sub.check_output(cmd,shell=True)).replace('\n',"")
        myimageen="data:image/png;base64,"+out
        return myimageen 
    def getImg(imgfld):
        myimage = glob.glob(imgfld+"/*.jpg")
        if len(myimage) == 0:
        	myimage = ""
        else:
            myimage = myimage[0]
        
        df = glob.glob(imgfld+"/*.tsv")
        	
        if len(df) == 0:
            return df, ""
      	    
        else:
            df = df[0]
        df = pd.read_csv(df, sep="\t")
        if len(list(df)) == 3:
    	    df = df.values
        else:
    	    df = [["No Genomic Islands Found","",""]]
        cgviewImg = imgProcess(myimage)
        return df, cgviewImg
        
    preResultDir = "/media/hd2/Data/genbank/all_JSCB/"
    imgfld = preResultDir+selected_value
    df, cgviewImg = getImg(imgfld)
    if len(df) == 0:
        return HttpResponse("Data Not Found. Please upload the genbank file for the genome and run the program again.")
    name = selected_value.replace("_"," ")
    mydict = {
        'cgviewImg':cgviewImg,
        'df':df,
        'name':name
    }
    template = loader.get_template (__package__+'/results.html')
    return HttpResponse(template.render(mydict,request))

def autocomplete_view(request):
    def are_all_substrings(substrings, target_string):
        for substring in substrings:
            if substring.lower() not in target_string.lower():
                return False
        return True
        
    query = request.GET['query'].split(" ")
    basedir = '/home/mehul/Mehul/djaenv/mysite/jscb/data/'
    f = open(basedir+'list_of_genomes')
    data = []
    for line in f:
        result = are_all_substrings(query, line.strip())
        if result: 
     	    data.append(line.strip()) 
    # Your logic to fetch and filter data for autocomplete
    #data = ["Option 1", "Option 2", "Option 3"]

    return JsonResponse(data, safe=False)
    
def main(request):
    def is_genbank_format(file_path):
        try:
            with open(file_path, 'r') as file:
                for line in file:
                    if line.startswith('LOCUS '):
                        return True
        except IOError:
            return False
        return False
        
    def imgProcess(imgpath):
        cmd = '''echo -n| base64 %s'''%(imgpath)
        out = sub.check_output(cmd,shell=True)
        print (type(out))
        myimageen="data:image/png;base64,"+str(out).replace('\n',"")
        """
        with open(imgpath,'rb') as img_file:
            myimageen = base64.b64encode(imgpath.encode('utf-8'))
            print (myimageen)
        """
        return myimageen 
        
    def getImg(imgfld):
        myimage = glob.glob(imgfld+"/*.jpg")[0]
        df = glob.glob(imgfld+"/*.tsv")[0]
        df = pd.read_csv(df, sep="\t")
        if len(list(df)) == 3:
    	    df = df.values
        else:
    	    df = [["No Genomic Islands Found","",""]]
        cgviewImg = imgProcess(myimage)
        return df, cgviewImg

    basedir = '/home/mehul/Mehul/djaenv/mysite/jscb/data/'
    testImage = imgProcess(basedir+"test")
    if request.method == 'POST':
        if request.POST.get("run"):
            cmd = '''bash /media/hd2/Data/jscb_docker/container_run.sh %s'''%(request.session["filename"])
            os.system(cmd)
            imgfld = basedir+request.session["filename"].replace(".","_")+"_JSCB"
            df, cgviewImg = getImg(imgfld)
    		
            mydict = {
                'cgviewImg':cgviewImg,
                'df':df
            }
            template = loader.get_template (__package__+'/results.html')
            return HttpResponse(template.render(mydict,request))

        form = FileUploadForm(request.POST, request.FILES)
        if form.is_valid():
            uploaded_file = form.cleaned_data['file']
            request.session["filename"] = uploaded_file.name
            with open(basedir + uploaded_file.name, 'wb+') as destination:
                for chunk in uploaded_file.chunks():
                    destination.write(chunk)
            if is_genbank_format(basedir + uploaded_file.name):
               print("The file is in GenBank format.")
            else:
                print("The file is not in GenBank format.")            
            if request.is_ajax():
                return JsonResponse({'message': 'File uploaded successfully'})
            else:
                return redirect('file_list')
    else:
        form = FileUploadForm()
        
    mydict = {
    	'testImage':testImage,
    	'form':form,
    }
    return render(request, __package__+'/file_upload.html', mydict)    

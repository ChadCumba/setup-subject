#!/usr/bin/env python
# encoding: utf-8
'''
preprocess.base -- shortdesc

preprocess.base is a description

It defines classes_and_methods

@author:     Chad Cumba 
            
@copyright:  2013 Imaging Research Center @ University of Texas. All rights reserved.
            
@license:    BSD

@contact:    chad.cumba@mail.utexas.edu
@deffield    updated: 22/10/13
'''


DEBUG = 1
TESTRUN = 0
PROFILE = 0
OUTPUT_DIR = '/work/01551/ccumba/output'
XNAT_SERVER = "https://xnat.irc.utexas.edu/xnat-irc"

import sys
if DEBUG:
    sys.path.insert(0, "/work/01551/ccumba/python_packages/nipype")
    sys.path.insert(0, "/work/01551/ccumba/setup_subject/src")
import os
import glob

"""
importing scipy or numpy causes extreme slowdown in the script, to avoid this so that the 
command line is snappy we only import modules that load numpy/scipy when they are needed.
listing them in this comment :

import xnat_tools

import nipype.pipeline.engine
import nipype.interfaces.utility
import nipype.interfaces.io
import nipype.interfaces.dcm2nii
import nipype.interfaces.fsl.preprocess
import nipype.interfaces.freesurfer.preprocess
import nipype.utils.filemanip
import nipype.interfaces.fsl.epi
from nipype.interfaces.freesurfer.preprocess import ReconAll

import qa
"""




from optparse import OptionParser

__all__ = []
__version__ = 0.1
__date__ = '2013-09-30'
__updated__ = '2013-09-30'




def name_it_bold(input):
    return 'bold'
    
def pop_last_get_dir(list_of_paths):
    import os
    if type(list_of_paths) is list:
        return os.path.split(list_of_paths[0])[0]
    return os.path.split(list_of_paths)[0]

def get_TR(dicom_header):
    import dicom
    import nipype.utils.filemanip
    
    dicom_header = nipype.utils.filemanip.loadpkl(dicom_header)
    
    return float(dicom_header.RepetitionTime/1000.0)
    

def get_dicom_headers(dicom_file):
    import dicom
    import nipype.utils.filemanip
    from nipype import logging
    
    if type(dicom_file) is list:
        dicom_file = dicom_file.pop()
    
    iflogger = logging.getLogger('interface')
    iflogger.debug('Getting headers from {}'.format(dicom_file))
    headers = dicom.read_file(dicom_file)
    
    
    
    iflogger.debug('Returning headers from {}'.format(dicom_file))
    nipype.utils.filemanip.savepkl(dicom_file + '.pklz', headers)
    
    return dicom_file + '.pklz'



def direct_nifti_to_directory(base_directory, dicom_header, niftis, bvals=None, bvecs=None):
    """Move nifti file to the correct directory for the subject
    @param dicom_headers: dicom_header object created by dicom.read_file and stored in a pickle dump
    @param nifti_file: a list of nifti files to move
    @param subject_directory: string representing the main directory for the subject

    @return nifti_destination: string representing where the file moved to

    """
    import dicom
    import os
    
    import nipype.utils.filemanip
    
    ANATOMY = 0
    BOLD = 1
    DTI = 2
    FIELDMAP = 3
    LOCALIZER = 4
    REFERENCE = 5
    DERIVED = 6
    
    from nipype import logging
    iflogger = logging.getLogger('interface')
    wlogger = logging.getLogger('workflow')
    
    iflogger.debug('Enetering direct nifti with inputs {} {} {} '.format(dicom_header, niftis, base_directory))
    
    dicom_header_filename = dicom_header
    dicom_header = nipype.utils.filemanip.loadpkl(dicom_header)
    
    scan_keys = [['MPRAGE','FSE','T1w','T2w','PDT','PD-T2','tse2d','t2spc','t2_spc'],
        ['epfid'],
        ['ep_b'],
        ['fieldmap','field_mapping'],
        ['localizer','Scout'],
        ['SBRef']]
        
    file_type = None 

    destination = ""
    
    if type(niftis) is str:
        niftis = [niftis]
    
    try:
        for i in range(0,6):
            for key in scan_keys[i]:
                if (dicom_header.ProtocolName.lower().find(key.lower()) > -1) or \
                    (dicom_header.SeriesDescription.lower().find(key.lower()) > -1) or \
                    (dicom_header.SequenceName.lower().find(key.lower()) > -1):
                        file_type = i
                        
        if dicom_header.ImageType[0] != "ORIGINAL":
            file_type = DERIVED
                
    except AttributeError, e:
        iflogger.warning("Nifti File(s) {} not processed because Dicom header dataset throwing error {} \n ".format(niftis, e)
                         + "Check your dicom headers")
        return []
    
        
    
    iflogger.debug('Direct nifti passed attribute checks')
    
    for nifti_file in niftis:
    
        iflogger.debug('classifying nifti_file {} from niftis {}'.format(nifti_file, niftis))
        nifti_basename = os.path.basename(nifti_file)
        
        nifti_extension = '.nii'
        
        if str(os.path.split(nifti_basename)[1]).endswith( '.gz'):
            nifti_extension = '.nii.gz'                  
        iflogger.debug('Passed filesplit on {}, file type is {}'.format(nifti_basename, file_type))
        
        try:
            run_number = nifti_basename.rsplit('a')[-2].rsplit('s')[-1].lstrip('0')
            run_name = dicom_header.ProtocolName.replace(' ','_')
               
        except Exception, e: 
            raise('unable to parse run number from nifti file {}'.format(nifti_file))

        if file_type == ANATOMY:
            try:
                run_number = nifti_basename.rsplit('a')[-2].rsplit('s')[-1].lstrip('0')
            except Exception, e:
                raise('unable to parse run number from nifti file {}'.format(nifti_file))
            
            mprage = False
            
            for key in ['mprage', 't1w']:
                if nifti_basename.lower().find(key) > 0:
                    mprage = True
            #MPRage files are high res anatomy and are named appropriately
            if nifti_basename.find('o') == 0 and mprage:
                
                highres_count = 0
                if os.path.exists(os.path.join(base_directory, 'anatomy')):
                    highres_count = len([ item for item in os.listdir(os.path.join(base_directory, 'anatomy')) if 'highres' in item])
                else:
                    os.makedirs(os.path.join(base_directory, 'anatomy'))
                destination = os.path.join(base_directory, 'anatomy', 'highres{0:03d}'.format(highres_count + 1)) + nifti_extension
                nipype.utils.filemanip.copyfile(
                    nifti_file, destination, copy=True)
                
            #PDT2 files are inplane brains
            if nifti_basename.lower().find('PDT2'.lower()) >= 0 and not mprage:
                inplane_count = 0
                if os.path.exists(os.path.join(base_directory, 'anatomy')):
                    inplane_count = len([ item for item in os.listdir(os.path.join(base_directory, 'anatomy')) if 'inplane' in item])
                else:
                    os.makedirs(os.path.join(base_directory, 'anatomy', 'other'))
                destination = os.path.join(base_directory, 'anatomy', 'inplane{0:03d}'.format(inplane_count + 1) + nifti_extension)
                nipype.utils.filemanip.copyfile(
                    nifti_file, destination, copy=True)
                nipype.utils.filemanip.copyfile(
                    nifti_file, os.path.join(base_directory, 'anatomy', 'other', nifti_basename), copy=True)
                
            else:
                #we aren't sure what type of anat file it is and just throw it in the "other" directory
                if not os.path.exists(os.path.join(base_directory, 'anatomy', 'other')):
                    os.makedirs(os.path.join(base_directory, 'anatomy', 'other'))
                destination = os.path.join(base_directory, 'anatomy', 'other', nifti_basename)
                nipype.utils.filemanip.copyfile(
                    nifti_file, destination, copy=True)
            
        elif file_type == BOLD:
            #bold files are stored in <dataset>/bold/runame_runnumber/bold.nii
            iflogger.debug("File {} is of type BOLD".format(nifti_file))
            try:
                run_number = nifti_basename.rsplit('a')[-2].rsplit('s')[-1].lstrip('0')
                run_name = dicom_header.ProtocolName.replace(' ','_')
                run_directory = os.path.join(base_directory, 'bold','%s_%s'%(run_name,run_number))
                
                if not os.path.exists(run_directory):
                    os.makedirs(run_directory)
                
            except Exception, e: 
                raise('unable to parse run number from nifti file {}'.format(nifti_file))
            
            destination = os.path.join(run_directory, 'bold' + nifti_extension)
            iflogger.debug("Sending file {} to destination {}".format(nifti_file, destination))
            nipype.utils.filemanip.copyfile(nifti_file, destination, copy=True)
            
        elif file_type == DTI:
            #DTI files are stored in <dataset>/dti/runname_runnnumber/DTI.nii.gz with a number starting at 1
            iflogger.debug('File type of {} is dti'.format(nifti_file))
            
            try:
                run_number = nifti_basename.rsplit('a')[-2].rsplit('s')[-1].lstrip('0')
                run_name = dicom_header.ProtocolName.replace(' ','_')
                run_directory = os.path.join(base_directory, 'dti','%s_%s'%(run_name,run_number))
                
                if not os.path.exists(os.path.join(base_directory, 'dti','%s_%s'%(run_name,run_number))):
                    os.makedirs(os.path.join(base_directory, 'dti','%s_%s'%(run_name,run_number)))
                
            except Exception, e: 
                raise('unable to parse run number from nifti file {}'.format(nifti_file))
            
            
            
            #copy the nifti into the dti dir with the rename
            destination = os.path.join(base_directory, 'dti', run_directory, 'dti' + nifti_extension)
            nipype.utils.filemanip.copyfile(nifti_file, destination, copy=True)
            
            #now for the bvecs and vals
            if type(bvecs) is not list:
                bvecs = [bvecs]
                
            for bvec in bvecs:
                destination = os.path.join(base_directory, 'dti', run_directory, "dti.bvec")
                nipype.utils.filemanip.copyfile(bvec, destination, copy=True)
            
            if type(bvals) is not list:
                bvals = [bvals]
            
            for bval in bvals:
                destination = os.path.join(base_directory, 'dti', run_directory, "dti.bval") 
                nipype.utils.filemanip.copyfile(bval, destination, copy=True)
            
            iflogger.debug('successfully moved DTI files {}, {}, {}'.format(nifti_file, bvecs, bvals))
        
        elif file_type == FIELDMAP:
            #here we're betting that the first fieldmap is the mag and the second is the phase
            #it's brittle but there isn't a good way to tell them apart from the headers.

            run_directory = os.path.join(base_directory, 'fieldmap', '{}_{}'.format(run_name, run_number))
            if not os.path.exists(run_directory):
                os.makedirs(run_directory)
            
            destination = os.path.join(run_directory, 'fieldmap' + nifti_extension)
            iflogger.debug('copying fieldmap {} to {}'.format(nifti_file, destination))
        
            nipype.utils.filemanip.copyfile(nifti_file, destination, copy=True)

        elif file_type == DERIVED:
            #derived filetypes we process the original data and leave it in the work directory
            #we don't really want or need them , but process them for 
            wlogger.info('File {} is a derived series, skipping further processing'.format(nifti_file))
            
    #this is a filetype we don't want to keep in openfmri style datasets, don't copy and return the original niftis
    #for compatibility
    if file_type == REFERENCE or file_type == LOCALIZER or file_type == DERIVED:
        return niftis

    #we didn't assign it a destination somehow even though it's a type of file we wanted to keep, throw an error
    if destination == "" :
        raise IOError('Failed to assign file {0} with filetype {1} to a directory'.format(niftis, file_type))
    else:
        #file is one we needed and was assigned a destination, copy its headers next to it
        root, ext = os.path.splitext(os.path.basename(destination))
        if root.endswith('nii'):
            root, ext = os.path.splitext(root)
        headers_destination = os.path.join(os.path.dirname(destination), root + '.pklz')
        nipype.utils.filemanip.copyfile(dicom_header_filename, headers_destination, copy=True)
        return niftis
    #this should never execute but will cause an error if we get here
    return destination
        
    

#@author Satrajit Ghosh
def get_subjectinfo(subject_id, base_dir, task_id, model_id):
    """Get info for a given subject

    Parameters
    ----------
    subject_id : string
        Subject identifier (e.g., sub001)
    base_dir : string
        Path to base directory of the dataset
    task_id : int
        Which task to process
    model_id : int
        Which model to process

    Returns
    -------
    run_ids : list of ints
        Run numbers
    conds : list of str
        Condition names
    TR : float
        Repetition time
    """
    from glob import glob
    import os
    import numpy as np
    condition_info = []
    cond_file = os.path.join(base_dir, 'models', 'model%03d' % model_id,
                             'condition_key.txt')
    with open(cond_file, 'rt') as fp:
        for line in fp:
            info = line.strip().split()
            condition_info.append([info[0], info[1], ' '.join(info[2:])])
    if len(condition_info) == 0:
        raise ValueError('No condition info found in %s' % cond_file)
    taskinfo = np.array(condition_info)
    n_tasks = len(np.unique(taskinfo[:, 0]))
    conds = []
    run_ids = []
    if task_id > n_tasks:
        raise ValueError('Task id %d does not exist' % task_id)
    for idx in range(n_tasks):
        taskidx = np.where(taskinfo[:, 0] == 'task%03d' % (idx + 1))
        conds.append([condition.replace(' ', '_') for condition
                      in taskinfo[taskidx[0], 2]])
        files = glob(os.path.join(base_dir,
                                  subject_id,
                                  'BOLD',
                                  'task%03d_run*' % (idx + 1)))
        run_ids.insert(idx, range(1, len(files) + 1))
    TR = np.genfromtxt(os.path.join(base_dir, 'scan_key.txt'))[1]
    return run_ids[task_id - 1], conds[task_id - 1], TR


def main(argv=None):
    '''Command line options.'''
    
    program_name = os.path.basename(sys.argv[0])
    program_version = "v0.1"
    program_build_date = "%s" % __updated__
 
    program_version_string = '%%prog %s (%s)' % (program_version, program_build_date)
    #program_usage = '''usage: spam two eggs''' # optional - will be autogenerated by optparse
    program_longdesc = '''''' # optional - give further explanation about what the program does
    program_license = "Copyright 2013 user_name (organization_name)                                            \
                Licensed under the Apache License 2.0\nhttp://www.apache.org/licenses/LICENSE-2.0"
 
    if argv is None:
        argv = sys.argv[1:]
    try:
        # setup option parser
        parser = OptionParser(version=program_version_string, epilog=program_longdesc, description=program_license)
        parser.add_option("-i", "--in", dest="infile", help="set input path [default: %default]", metavar="FILE")
        parser.add_option("-o", "--out", dest="outfile", help="set output path [default: %default]", metavar="FILE")
        parser.add_option("-v", "--verbose", dest="verbose", action="count", help="set verbosity level [default: %default]")
        parser.add_option("-d", "--datadir", dest="data_directory", help="set project data directory [default: %default/ ]",
                          metavar="DIRECTORY")
        parser.add_option("-s", "--subject", dest="subject", help="set subject directory [default: %default]", 
                          metavar="DIRECTORY")
        parser.add_option("-m", "--model", dest="model_id", help="set the model to perform analysis on [default: %default]",
                          metavar="NUMBER")
        parser.add_option("-t", "--task", dest="task_id", help="set the task to perform analysis on [default: %default]",
                          metavar="NUMBER")
        parser.add_option("-w", "--work", dest="work_directory", help="set the directory to store transitional files [default: %default]",
                          metavar="DIRECTORY")
        parser.add_option("-g", "--getdata", dest="get_data", help="get data from XNAT [default: %default]", action="store_true")
        parser.add_option("-p", "--project", dest="project", help="set the XNAT project name [default: %default]", 
                          metavar="XNAT PROJECT")
        parser.add_option("--dcm2nii", dest="dicom_to_nifti", help="convert dicom files to nifti format", action="store_true")
        parser.add_option("--motcorr", dest="motion_correction", help="run bet motion correction", action="store_true")
        parser.add_option("--melodic", dest="melodic", help="run melodic on func data", action="store_true")
        parser.add_option("--betfunc", dest="skull_strip", help="run BET on func data", action="store_true")
        parser.add_option("--fsrecon", dest="autorecon_all", help="run freesurfer autorecon1", action="store_true")
        parser.add_option("--dtiqa", dest="dti_qa", help="generate QA report on DTI data", action="store_true")
        parser.add_option("--qa", dest="fmri_qa", help="generate QA report on bold data", action="store_true")
        parser.add_option("--fm", dest="process_fieldmap", help="process fieldmap files", action="store_true")
        
         
        # set defaults
        parser.set_defaults(outfile="./out.txt", infile="./in.txt", work_directory=os.environ.get('SCRATCH'), task_id=1,
                            model_id=1, subject="all", get_data=False, data_directory='.', project=None, motion_correction=False,
                            melodic=False, skull_strip=False, autorecon_all=False, dti_qa=False, fmri_qa=False,
                            process_fieldmap=False)
        # process options
        (opts, args) = parser.parse_args(argv)
        if opts.verbose > 0:
            print("verbosity level = %d" % opts.verbose)
        if opts.infile:
            print("infile = %s" % opts.infile)
        if opts.outfile:
            print("outfile = %s" % opts.outfile)
        if opts.subject == "all":
            opts.subject = None
        if opts.get_data:
            if opts.project is None:
                raise IOError('Must specify --project when using --getdata')
            if opts.subject == "all" or opts.subject is None:
                raise IOError('Must specify --subject when using --getdata, and --subject cannot be "all"')
        opts.data_directory = os.path.abspath(opts.data_directory)
        
        # MAIN BODY #
    except Exception, e:
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help\n")
        return 2
    
    if opts.get_data :
        import xnat_tools
        #we aren't doing this with an XNATSource node because the current implementation cannot download
        #more than 10 files at once
        if DEBUG:
            print XNAT_SERVER
            print os.path.join(opts.data_directory, 'raw')
            print opts.project
            print opts.subject
        print "Type your XNAT username and password below"
        if not os.path.exists(os.path.join(opts.data_directory, opts.subject, 'raw')):
                os.makedirs(os.path.join(opts.data_directory, opts.subject, 'raw'))
        xnat_tools.down_subject_dicoms(XNAT_SERVER, os.path.join(opts.data_directory, opts.subject, 'raw'), opts.project, opts.subject)

    

    if opts.dicom_to_nifti:
        if opts.subject == "all":
            for subject in os.listdir(opts.data_directory):
                workflow = openfmri_dicom_to_nifti(opts.data_directory, subject)
                workflow.run()
        else:
            workflow = openfmri_dicom_to_nifti(opts.data_directory, opts.subject)
            workflow.run()
    if any([opts.autorecon_all, opts.melodic, opts.motion_correction, opts.skull_strip, opts.dti_qa, opts.fmri_qa, opts.process_fieldmap]):
        workflow = preprocess_dataset(opts.data_directory, opts.subject, model_id=opts.model_id, task_id=opts.task_id,
                                 work_directory=opts.work_directory, xnat_datasource=opts.get_data,
                                 run_motion_correction=opts.motion_correction, run_skull_strip=opts.skull_strip, 
                                 run_melodic=opts.melodic, run_autorecon=opts.autorecon_all, run_dti_qa=opts.dti_qa,
                                 run_fmri_qa=opts.fmri_qa, run_process_fieldmap=opts.process_fieldmap)
        workflow.run()

def preprocess_dataset(data_directory, subject, model_id=None, task_id=None,  work_directory=None, xnat_datasource=False,
                             run_motion_correction=False, run_skull_strip=False, run_melodic=False, run_autorecon=False,
                             run_dti_qa=False, run_fmri_qa=False, run_process_fieldmap=False):
    """
    @TODO - docs
    
    @param data_directory - directory containing subject folders
    @param subject - directory within data_directory containing a single subjects data
    @param model_id - the model number to run
    @param task_id - the task to run
    
    @param work_directory - directory to store transitional files
    """
    
    import nipype.pipeline.engine
    import nipype.interfaces.utility
    import nipype.interfaces.io
    import nipype.interfaces.fsl.preprocess
    import nipype.interfaces.fsl.epi
    
    import qa
    
    from nipype.interfaces.freesurfer.preprocess import ReconAll 

    from nipype import logging
    wflogger = logging.getLogger('workflow')

    if not os.path.exists(os.path.join(data_directory, subject)):
        raise IOError("Subject directory {0} not found".format(os.path.join(data_directory,subject)))
    
    subject_directory = os.path.join(data_directory, subject)
    
    bold_directories = [path.split(os.path.sep)[-1] for path in glob.glob(os.path.join(subject_directory, '[Bb][Oo][Ll][Dd]', '*'))]
    dti_directories = [path.split(os.path.sep)[-1] for path in glob.glob(os.path.join(subject_directory, '[Dd][Tt][Ii]', '*'))]
    fieldmap_directories = [path.split(os.path.sep)[-1] for path in glob.glob(os.path.join(subject_directory, '[Ff]ieldmap', '*'))]
        
    wflogger.debug('prepping infosource')
    infosource = nipype.pipeline.engine.Node(
                    nipype.interfaces.utility.IdentityInterface(fields=["bold_dirs", "dti_dirs", "fieldmap_dirs"]),
                    name="infosource")
    
    
    #setup the datasources
    infosource.iterables = []
    infosource.iterables.append(("bold_dirs", bold_directories))
    boldsource = nipype.pipeline.engine.Node(
                        nipype.interfaces.io.DataGrabber(infields=["bold_dirs"],
                                                         outfields=["bold", "headers"]),
                                            name="boldsource")
    boldsource.inputs.base_directory = os.path.join(subject_directory, "bold")
    boldsource.inputs.template = "*"
    boldsource.inputs.sort_filelist = True
    boldsource.inputs.field_template = {"bold" : "%s/bold.nii*"}
    boldsource.inputs.template_args = {"bold" : [["bold_dirs"]]}
    
    workflow = nipype.pipeline.engine.Workflow(name="openfmri")
    
    if work_directory is None:
        work_directory = os.path.join(os.getcwd(),'working')
        
    workflow.base_dir = work_directory
    workflow.connect(infosource, "bold_dirs", boldsource, "bold_dirs")
    #connect iterable directories if they exist
    
    #Data storage and sink
    datasink = nipype.pipeline.engine.Node(nipype.interfaces.io.DataSink(), name="datasink")
    
    datasink.inputs.base_directory = subject_directory
    
    if DEBUG:
        datasink.inputs.base_directory = OUTPUT_DIR
    
    #FSL mcflirt motion correction
    if any([run_motion_correction, run_skull_strip, run_melodic, run_fmri_qa]):
        motion_correction = nipype.pipeline.engine.Node(
                                interface=nipype.interfaces.fsl.preprocess.MCFLIRT(),
                                name="motion_correction")
        motion_correction.inputs.terminal_output = "file"
        workflow.connect(boldsource, "bold", motion_correction, "in_file")
        workflow.connect(motion_correction, "out_file", datasink, 'bold.@bold_mcf')

    #run QA report for dti files
    if run_dti_qa:
        
        
        dti_qa = nipype.pipeline.engine.Node(
                    interface=qa.DTIQATask(),
                    name="dti_qa")
        
        workflow.connect(dti_qa, "snr", datasink, "dti_qa")
        workflow.connect(dti_qa, "slice_corr", datasink, "dti_qa.@slice_corr")
        workflow.connect(dti_qa, "interleave_corr", datasink, "dti_qa.@interleave_corr")
        workflow.connect(dti_qa, "fa", datasink, "dti_qa.@fa")
        workflow.connect(dti_qa, "fd", datasink, "dti_qa.@fd")
        workflow.connect(dti_qa, "sse", datasink, "dti_qa.@sse")
        workflow.connect(dti_qa, "worst_grad", datasink, "dti_qa.@worst_grad")
        workflow.connect(dti_qa, "report", datasink, "dti_qa.@report")
    
    if run_fmri_qa:
        
        boldsource.inputs.field_template['headers'] = "%s/bold.pklz"
        boldsource.inputs.template_args['headers'] = [[ "bold_dirs"]]
        
        temporal_resolution = nipype.pipeline.engine.Node(
                                nipype.interfaces.utility.Function(input_names=["dicom_header"],
                                                                   output_names=["temporal_resolution"],
                                                                   function=get_TR),
                                                          name="temporal_resolution")
        
        fmri_qa = nipype.pipeline.engine.Node(
                    interface=qa.FMRIQATask(),
                    name="fmri_qa")
        
        workflow.connect(motion_correction, "out_file", fmri_qa, "in_file")
        workflow.connect(boldsource, "headers", temporal_resolution, "dicom_header")
        workflow.connect(temporal_resolution, "temporal_resolution", fmri_qa, "temporal_resolution")
        
        workflow.connect(fmri_qa, "confound", datasink, "fmri_qa.@confound") 
        workflow.connect(fmri_qa, "dvars_png", datasink, "fmri_qa.@dvars_png")
        workflow.connect(fmri_qa, "dvars_list", datasink, "fmri_qa.@dvars_list")
        workflow.connect(fmri_qa, "fd_list" ,datasink, "fmri_qa.@fd_list")
        workflow.connect(fmri_qa, "fd_png" ,datasink, "fmri_qa.@fd_png")
        workflow.connect(fmri_qa, "mad" ,datasink, "fmri_qa.@mad")
        workflow.connect(fmri_qa, "maskmean",datasink, "fmri_qa.@maskmean")
        workflow.connect(fmri_qa, "meanpsd" ,datasink, "fmri_qa.@meanpsd")
        workflow.connect(fmri_qa, "qadata" ,datasink, "fmri_qa.@qadata")
        workflow.connect(fmri_qa, "qa_report",datasink, "fmri_qa.@qa_report")
        workflow.connect(fmri_qa, "spike" ,datasink, "fmri_qa.@spike")
        workflow.connect(fmri_qa, "voxcv" ,datasink, "fmri_qa.@voxcv")
        workflow.connect(fmri_qa, "voxmean",datasink, "fmri_qa.@voxmean")
        workflow.connect(fmri_qa, "voxsfnr_volume",datasink, "fmri_qa.@voxsfnr_volume")
        workflow.connect(fmri_qa, "voxsfnr_png" ,datasink, "fmri_qa.@voxsfnr_png")
    
    
    #betfunc skull stripping
    if run_skull_strip:
        skull_strip = nipype.pipeline.engine.Node(
                        interface=nipype.interfaces.fsl.BET(),
                        name="skull_strip")
        workflow.connect(motion_correction, "out_file", skull_strip, "in_file")
        skull_strip.inputs.terminal_output = 'file'
        skull_strip.inputs.functional = True
        workflow.connect(skull_strip, "out_file", datasink, "bold.@bold_mcf_brain")
        
    #melodic analysis
    if run_melodic:
        independent_components_analysis = nipype.pipeline.engine.Node(
                                            interface=nipype.interfaces.fsl.model.MELODIC(),
                                            name="melodic")
        workflow.connect(motion_correction,"out_file", independent_components_analysis, "in_files")
        independent_components_analysis.inputs.terminal_output = "file"
        workflow.connect(independent_components_analysis, "out_dir", datasink, "melodic")
        workflow.connect(independent_components_analysis, "report_dir", datasink, "melodic_report")
    
    #freesurfer autorecon all
    if run_autorecon:
        cortical_reconstruction = nipype.pipeline.engine.Node(
                                    interface=ReconAll(),
                                    name="cortical_reconstruction")
        
        workflow.connect(cortical_reconstruction, 'subject_id', datasink, 'container')
        workflow.connect(cortical_reconstruction, 'T1', datasink, 'highres')
    
    #fieldmap preprocessing
    if run_process_fieldmap:
        
        prepare_fieldmap = nipype.pipeline.engine.Node(
                            interface=nipype.interfaces.fsl.epi.PrepareFieldmap(),
                            name="prepare_fieldmap")
        prepare_fieldmap.inputs.output_type="NIFTI_GZ"
        
        
        
        #skull strip the fieldmap file
        fieldmap_skull_strip = nipype.pipeline.engine.Node(
                        interface=nipype.interfaces.fsl.BET(),
                        name="fieldmap_skull_strip")
        fieldmap_skull_strip.inputs.terminal_output = 'file'
                
        workflow.connect(fieldmap_skull_strip, "out_file", prepare_fieldmap, "in_phase")
        workflow.connect(fieldmap_skull_strip, "out_file", prepare_fieldmap, "in_magnitude")
        workflow.connect(prepare_fieldmap, "out_fieldmap", datasink, "fieldmap")
    
    
    return workflow

def analyze_openfmri_dataset(data_directory, subject=None, model_id=None, task_id=None,  work_directory=None, xnat_datasource=False,
                             run_motion_correction=False, run_skull_strip=False, run_melodic=False, run_autorecon=False, 
                             run_bet_inplane=False):
    """
    @TODO - docs
    
    @param data_directory - directory containing subject folders
    @param subject - directory within data_directory containing a single subjects data
    @param model_id - the model number to run
    @param task_id - the task to run
    
    @param work_directory - directory to store transitional files
    """
    subjects = [path.split(os.path.sep)[-1] for path in glob.glob(os.path.join(data_directory, '*'))]
    
    infosource = nipype.pipeline.engine.Node(
                    nipype.interfaces.utility.IdentityInterface(fields=["subject_id",
                                                                        "model_id",
                                                                        "task_id"]),
                                                                name="infosource")
    
    if subject is None:
        infosource.iterables = [("subject_id", subjects),
                                ("model_id", [model_id])]
    else:
        infosource.iterables = [("subject_id", [subjects[subjects.index(subject)]]),
                                ("model_id", [model_id]),
                                ("task_id", [task_id]),
                                ]
    
    subject_info = nipype.pipeline.engine.Node(
                        nipype.interfaces.utility.Function(input_names=["subject_id", "base_dir",
                                                                        "task_id", "model_id"],
                                                           output_names=["run_id", "conds", "TR"],
                                                           function=get_subjectinfo),
                                               name="subjectinfo")
    
    subject_info.inputs.base_dir = data_directory

    #setup the datasources
    datasource = nipype.pipeline.engine.Node(
                    nipype.interfaces.io.DataGrabber(infields=["subject_id", "run_id",
                                                               "task_id", "model_id"],
                                                     outfields=["anat", "bold", "behav"]),
                                             name="datasource")
    
    datasource.inputs.base_directory = data_directory
    datasource.inputs.template = "*"
    datasource.inputs.field_template = {"anat" : "%s/anatomy/highres001.nii.gz",
                                        "bold" : "%s/BOLD/task%03d_run%03d/bold.nii.gz",
                                        "behav" : "%s/model/model%03d/onsets/task%03d_run%03d/cond*.txt"
                                        }
    datasource.inputs.template_args = {"anat" : [["subject_id"]],
                                       "bold" : [["subject_id", "task_id", "run_id"]],
                                       "behav" : [["subject_id", "model_id", "task_id", "run_id"]]
                                       }
    datasource.inputs.sort_filelist = True
    #This turns off the errors if a grabber cannot find a file
    #Note that the nodes later on will raise if they are executed and do not
    #have proper inputs, so it is not necessary
    datasource.inputs.raise_on_empty = False
    
    workflow = nipype.pipeline.engine.Workflow(name="openfmri")
    
    if work_directory is None:
        work_directory = os.path.join(os.getcwd(),'working')
        
    workflow.base_dir = work_directory
    
    workflow.connect(infosource, "subject_id", subject_info, "subject_id")
    workflow.connect(infosource, "model_id", subject_info, "model_id")
    workflow.connect(infosource, "task_id", subject_info, "task_id")
    workflow.connect(infosource, "subject_id", datasource, "subject_id")
    workflow.connect(infosource, "model_id", datasource, "model_id")
    workflow.connect(infosource, "task_id", datasource, "task_id")    
    workflow.connect(subject_info, 'run_id', datasource, 'run_id')
    #Data storage and sink
    datasink = nipype.pipeline.engine.Node(nipype.interfaces.io.DataSink(), name="datasink")
    datasink.inputs.base_directory = OUTPUT_DIR
    
    #FSL mcflirt motion correction
    if any([run_motion_correction, run_skull_strip, run_melodic]):
        motion_correction = nipype.pipeline.engine.Node(
                                interface=nipype.interfaces.fsl.preprocess.MCFLIRT(),
                                name="motion_correction")
        motion_correction.inputs.terminal_output = "file"
        workflow.connect(datasource, "bold", motion_correction, "in_file")
        workflow.connect(motion_correction, "out_file", datasink, 'bold_mcf')

    #run bet on inplane brain
    if run_bet_inplane:
       inplane_skull_strip = nipype.pipeline.engine.Node(
                                interface=nipype.interfaces.fsl.BET(),
                                name="inplane_skull_strip")
       inplane_skull_strip.inputs.terminal_output = "file"
       workflow.connect(datagrabber, "inplane", skull_strip, "in_file")
       workflow.connect(skull_strip, "out_file", datasink, "inplane_brain")
    
    #betfunc skull stripping
    if run_skull_strip:
        skull_strip = nipype.pipeline.engine.Node(
                        interface=nipype.interfaces.fsl.BET(),
                        name="skull_strip")
        workflow.connect(motion_correction, "out_file", skull_strip, "in_file")
        skull_strip.inputs.terminal_output = 'file'
        workflow.connect(skull_strip, "out_file", datasink, "bold_mcf_brain")
        
    #melodic analysis
    if run_melodic:
        independent_components_analysis = nipype.pipeline.engine.Node(
                                            interface=nipype.interfaces.fsl.model.MELODIC(),
                                            name="melodic")
        workflow.connect(motion_correction,"out_file", independent_components_analysis, "in_files")
        independent_components_analysis.inputs.terminal_output = "file"
        workflow.connect(independent_components_analysis, "out_dir", datasink, "melodic")
        workflow.connect(independent_components_analysis, "report_dir", datasink, "melodic_report")
    
    #freesurfer autorecon all
    if run_autorecon:
        cortical_reconstruction = nipype.pipeline.engine.Node(
                                    interface=ReconAll(),
                                    name="cortical_reconstruction")
        
        workflow.connect(datasource, "anat", cortical_reconstruction, "T1_files" )
        workflow.connect(cortical_reconstruction, 'subject_id', datasink, 'container')
        workflow.connect(cortical_reconstruction, 'T1', datasink, 'highres')
    

    
    
    
    
    workflow.run(plugin="SGE", plugin_args={"qsub_args":("-l h_rt=3:00:00 -q normal -A Analysis_Lonestar " 
                                            "-pe 12way 24 -M chad.cumba@mail.utexas.edu")})

def openfmri_dicom_to_nifti(openfmri_subject_directory, subject_id):
    
    import nipype.pipeline.engine
    import nipype.interfaces.utility
    import nipype.interfaces.io
    import nipype.interfaces.dcm2nii
    
    scans_container = os.path.join(openfmri_subject_directory, subject_id, 'raw', subject_id)
    scan_directories = [directory for directory in os.listdir(scans_container)  ]
    #this is all just to tell it to iterate over the scan directories
    infosource = nipype.pipeline.engine.Node(
                    interface=nipype.interfaces.utility.IdentityInterface(fields=['scan_id']),
                    name="infosource")

    infosource.iterables = ('scan_id', scan_directories)
    
    #here we're telling the rest of the nodes where things exist on disk
    datasource = nipype.pipeline.engine.Node(
                    nipype.interfaces.io.DataGrabber(infields=["scan_id"],
                                                     outfields=["dicoms"]),
                                             name="datasource")
    datasource.inputs.base_directory = scans_container
    datasource.inputs.template = "*dcm"
    datasource.inputs.field_template = {'dicoms' : '%s/*dcm'}
    datasource.inputs.template_args = {'dicoms' : [['scan_id']]}
    datasource.inputs.sorted = False
    datasource.inputs.sort_filelist = False
    datasource.inputs.ignore_exception = True

    
    #build the nifti converter
    converter = nipype.pipeline.engine.Node(
                    nipype.interfaces.dcm2nii.Dcm2nii(terminal_output="file"),
                    name="nifticonvert"
    )

    #force 4d volumes and no dates in filenames
    converter.inputs.args = '-4 y -d n'
    converter.inputs.nii_output = True
    #ignore exceptions on bad inputs to prevent a full crash from one bad scan
    converter.inputs.ignore_exception = True
    converter.inputs.gzip_output = True
    
    #this is a function interface that grabs the dicom info for later
    dicomheaders = nipype.pipeline.engine.Node(
                        nipype.interfaces.utility.Function(input_names=["dicom_file"],
                                                           output_names=["dicom_header"],
                                                           function=get_dicom_headers),
                                               name="dicomheaders")


    openfmri_organization = nipype.pipeline.engine.Node(
                                nipype.interfaces.utility.Function(input_names=["base_directory", "dicom_header","niftis","bvecs", "bvals" ],
                                                                   output_names = ["destination"],
                                                                   function=direct_nifti_to_directory),
                                                        name="openfmri_organization")
    openfmri_organization.inputs.base_directory = os.path.join(openfmri_subject_directory, subject_id)

    #just a test sink for now
    datasink = nipype.pipeline.engine.Node(
                    nipype.interfaces.io.DataSink(base_directory=os.path.join(openfmri_subject_directory, subject_id)),
                    name="datasink")

    datasink.inputs.substitutions = [('_scan_id_', '')] 
    datasink.inputs.ignore_exception = True
    
    workflow = nipype.pipeline.engine.Workflow(name='dicom2nifti')
    workflow.connect(infosource, 'scan_id', datasource, 'scan_id')
    #lambda function that just runs pop() on any list you hand it
    pop_last = lambda x: [x[0]] if type(x) is list else x
    #workflow.connect(datasource, ('dicoms', pop_last_get_dir), converter, "dicom_dir")
    workflow.connect(datasource, ('dicoms', pop_last), converter, "source_names")
    workflow.connect(datasource, ('dicoms', pop_last), dicomheaders, "dicom_file")
    workflow.connect(dicomheaders, "dicom_header", openfmri_organization, "dicom_header")
    workflow.connect(converter, "converted_files", openfmri_organization, "niftis")
    workflow.connect(converter, "bvals", openfmri_organization, "bvals")
    workflow.connect(converter, "bvecs", openfmri_organization, "bvecs")
    workflow.connect(converter, "converted_files", datasink, "niftis")

    workflow.base_dir = os.environ.get("SCRATCH", "/tmp")
    return workflow
    

if __name__ == "__main__":
    if DEBUG:
        pass
    if TESTRUN:
        import doctest
        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'preprocess.base_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())


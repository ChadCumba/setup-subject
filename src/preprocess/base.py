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

import sys
import os
import glob
import xnat_tools
import dicom

import nipype.pipeline.engine
import nipype.interfaces.utility
import nipype.interfaces.io
import nipype.interfaces.dcm2nii
import nipype.interfaces.fsl.preprocess
import nipype.utils.filemanip


from nipype.interfaces.freesurfer.preprocess import ReconAll

from optparse import OptionParser

__all__ = []
__version__ = 0.1
__date__ = '2013-09-30'
__updated__ = '2013-09-30'

DEBUG = 1
TESTRUN = 0
PROFILE = 0
OUTPUT_DIR = '/work/01551/ccumba/output'
XNAT_SERVER = "https://xnat.irc.utexas.edu/xnat-irc"



def name_it_bold(input):
    return 'bold'
    
def pop_last_get_dir(list_of_paths):
    import os
    if type(list_of_paths) is list:
        return os.path.split(list_of_paths.pop())[0]
    return os.path.split(list_of_paths)[0]

def get_dicom_headers(dicom_file):
    import dicom
    import pickle
    from nipype import logging
    iflogger = logging.getLogger('interface')
    iflogger.debug('Getting headers from {}'.format(dicom_file))
    headers = dicom.read_file(dicom_file)
    
    
    
    iflogger.debug('Returning headers from {}'.format(dicom_file))
    
    with open(dicom_file + '.pklz', 'w') as pickle_handle:
        pickle.dump(headers, pickle_handle)
    
    return dicom_file + '.pklz'

def direct_nifti_to_directory(dicom_header, niftis, base_directory):
    """Move nifti file to the correct directory for the subject
    @param dicom_headers: dicom_header object created by dicom.read_file and stored in a pickle dump
    @param nifti_file: a list of nifti files to move
    @param subject_directory: string representing the main directory for the subject

    @return nifti_destination: string representing where the file moved to

    """
    import dicom
    import pickle
    import os
    
    import nipype.utils.filemanip
    
    ANATOMY = 0
    BOLD = 1
    DTI = 2
    FIELDMAP = 3
    LOCALIZER = 4
    REFERENCE = 5
    
    from nipype import logging
    iflogger = logging.getLogger('interface')
    
    iflogger.debug('Enetering direct nifti with inputs {} {} {} '.format(dicom_header, niftis, base_directory))
    
    with open(dicom_header, 'r') as dicom_headers_handle:
        dicom_header = pickle.load(dicom_headers_handle)
    
    scan_keys = [['MPRAGE','FSE','T1w','T2w','PDT','PD-T2','tse2d','t2spc','t2_spc'],
        ['epfid'],
        ['ep_b'],
        ['fieldmap','field_mapping'],
        ['localizer','Scout'],
        ['SBRef']]
        
    file_type = None 

    destination = []
    
    if type(niftis) is str:
        niftis = [niftis]
    
    try:
        for i in range(0,6):
            for key in scan_keys[i]:
                if (dicom_header.ProtocolName.lower().find(key.lower()) > -1) or \
                    (dicom_header.SeriesDescription.lower().find(key.lower()) > -1) or \
                    (dicom_header.SequenceName.lower().find(key.lower()) > -1):
                        file_type = i
    except AttributeError, e:
        iflogger.warning("Nifti File(s) {} not processed because Dicom header dataset throwing error {} \n ".format(niftis, e)
                         + "Check your dicom headers")
        return []
    
    iflogger.debug('Direct nifti passed attribute checks')
    
    for nifti_file in niftis:
    
        iflogger.debug('classifying nifti_file {} from niftis {}'.format(nifti_file, niftis))
        nifti_basename = os.path.basename(nifti_file)
        
        nifti_extension = '.nii'
        
        if os.path.split(nifti_basename)[1] == '.gz':
            nifti_extension = '.nii.gz'                  
        iflogger.debug('Passed filesplit on {}, file type is {}'.format(nifti_basename, file_type))
        
        if file_type == ANATOMY:
            try:
                run_number = nifti_basename.rsplit('a')[-2].rsplit('s')[-1].lstrip('0')
            except Exception, e:
                raise('unable to parse run number from nifti file {}'.format(nifti_file))
            
            mprage = False
            
            for key in ['mprage', 't1w']:
                if nifti_basename.lower().find(key) > 0:
                    mprage = True
            
            if nifti_basename.find('o') == 0 and mprage:
                
                highres_count = 0
                if os.path.exists(os.path.join(base_directory, 'anatomy')):
                    highres_count = len([ item for item in os.listdir(os.path.join(base_directory, 'anatomy')) if 'highres' in item])
                else:
                    os.makedirs(os.path.join(base_directory, 'anatomy'))
                destination.append(os.path.join(base_directory, 'anatomy', 'highres{0:03d}'.format(highres_count + 1)) + nifti_extension)
                nipype.utils.filemanip.copyfile(
                    nifti_file, destination[-1], copy=True)
                
            if nifti_basename.lower().find('PDT2'.lower()) > 0 and not mprage:
                inplane_count = 0
                if os.path.exists(os.path.join(base_directory, 'anatomy')):
                    inplane_count = len([ item for item in os.listdir(os.path.join(base_directory, 'anatomy')) if 'inplane' in item])
                else:
                    os.makedirs(os.path.join(base_directory, 'anatomy', 'other'))
                destination.append(os.path.join(base_directory, 'anatomy', 'inplane{0:03d}'.format(inplane_count + 1) + nifti_extension))
                nipype.utils.filemanip.copyfile(
                    nifti_file, destination[-1], copy=True)
                nipype.utils.filemanip.copyfile(
                    nifti_file, os.path.join(base_directory, 'anatomy', 'other', nifti_basename), copy=True)
                
            else:
                if not os.path.exists(os.path.join(base_directory, 'anatomy', 'other')):
                    os.makedirs(os.path.join(base_directory, 'anatomy', 'other'))
                destination.append(os.path.join(base_directory, 'anatomy', 'other', nifti_basename))
                nipype.utils.filemanip.copyfile(
                    nifti_file, destination[-1], copy=True)
            
        elif file_type == BOLD:
            iflogger.debug("File {} is of type BOLD".format(nifti_file))
            try:
                run_number = nifti_basename.rsplit('a')[-2].rsplit('s')[-1].lstrip('0')
                run_name = dicom_header.ProtocolName.replace(' ','_')
                run_directory = os.path.join(base_directory, 'bold','%s_%s'%(run_name,run_number))
                
                if not os.path.exists(os.path.join(base_directory, 'bold','%s_%s'%(run_name,run_number))):
                    os.makedirs(os.path.join(base_directory, 'bold','%s_%s'%(run_name,run_number)))
                
            except Exception, e: 
                raise('unable to parse run number from nifti file {}'.format(nifti_file))
            
            destination.append(os.path.join(base_directory, 'bold', run_directory, 'bold' + nifti_extension))
            iflogger.debug("Sending file {} to destination {}".format(nifti_file, destination))
            nipype.utils.filemanip.copyfile(nifti_file, destination[-1], copy=True)
            
        elif file_type == DTI:
            
            iflogger.debug('File type of {} is dti'.format(nifti_file))
            
            dti_count = 0
            if os.path.exists(os.path.join(base_directory, 'dti')):
                iflogger.debug('dit path {} exists, skipping creation'.format(os.path.join(base_directory, 'dti')))
                dti_count = int(len([ item for item in os.listdir(os.path.join(base_directory, 'dti')) if 'DTI' in item ]) / 3)
            else:
                iflogger.debug('Creating dti directory {}'.format(os.path.join(base_directory, 'dti')))
                os.makedirs(os.path.join(base_directory, 'dti'))
            
            destination.append(os.path.join(base_directory, 'dti', 'DTI_{0:03d}'.format(dti_count + 1) + nifti_extension))
            nipype.utils.filemanip.copyfile(nifti_file, destination[-1], copy=True)
        
        elif file_type == FIELDMAP:
            fieldmap_count = 0
            if os.path.exists(os.path.join(base_directory, 'fieldmap')):
                fieldmap_count = len([ item for item in os.listdir(os.path.join(base_directory, 'fieldmap')) if 'fieldmap' in item ])
                iflogger.debug('fieldmap count is {}'.format(fieldmap_count))
            else:
                os.makedirs(os.path.join(base_directory, 'fieldmap'))
        
            if fieldmap_count == 0:
                destination.append(os.path.join(base_directory, 'fieldmap', 'fieldmap_mag' + nifti_extension))
                nipype.utils.filemanip.copyfile(nifti_file, destination[-1], copy=True)
            elif fieldmap_count == 1:
                destination.append(os.path.join(base_directory, 'fieldmap', 'fieldmap_phase' + nifti_extension))
                nipype.utils.filemanip.copyfile(nifti_file, destination[-1], copy=True)
            else:
                destination.append(os.path.join(base_directory, 'fieldmap', 'fieldmap_{0:03d}'.format(fieldmap_count + 1) + nifti_extension))
                nipype.utils.filemanip.copyfile(nifti_file, destination[-1], copy=True)

    if destination == [] and not (file_type == REFERENCE or file_type == LOCALIZER):
        raise IOError('Failed to assign file {0} with filetype {1} to a directory'.format(niftis, file_type))
    else:
        return niftis
    
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
         
        # set defaults
        parser.set_defaults(outfile="./out.txt", infile="./in.txt", work_directory=os.environ.get('SCRATCH'), task_id=1,
                            model_id=1, subject="all", get_data=False, data_directory='.', project=None, motion_correction=False,
                            melodic=False, skull_strip=False, autorecon_all=False)
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
            [openfmri_dicom_to_nifti(opts.data_directory, subject) for subject in os.listdir(opts.data_directory)]
        else:
            openfmri_dicom_to_nifti(opts.data_directory, opts.subject)
    if any([opts.autorecon_all, opts.melodic, opts.motion_correction, opts.skull_strip]):
        preprocess_dataset(opts.data_directory, subject=opts.subject, model_id=opts.model_id, task_id=opts.task_id,
                                 work_directory=opts.work_directory, xnat_datasource=opts.get_data,
                                 run_motion_correction=opts.motion_correction, run_skull_strip=opts.skull_strip, 
                                 run_melodic=opts.melodic, run_autorecon=opts.autorecon_all)

def preprocess_dataset(data_directory, subject=None, model_id=None, task_id=None,  work_directory=None, xnat_datasource=False,
                             run_motion_correction=False, run_skull_strip=False, run_melodic=False, run_autorecon=False):
    """
    @TODO - docs
    
    @param data_directory - directory containing subject folders
    @param subject - directory within data_directory containing a single subjects data
    @param model_id - the model number to run
    @param task_id - the task to run
    
    @param work_directory - directory to store transitional files
    """
    subjects = [path.split(os.path.sep)[-1] for path in glob.glob(os.path.join(data_directory, '*'))]
    
    bolddirs_by_subject = {}
    for item in subjects:
        bolddirs_by_subject[item] = [path.split(os.path.sep)[-1] for path in glob.glob(os.path.join(data_directory, item, 'bold', '*'))]
    
    
    infosource = nipype.pipeline.engine.Node(
                    nipype.interfaces.utility.IdentityInterface(fields=["subject_id", "bold_dirs"]),
                                                                name="infosource")
    
    if subject is None:
        infosource.iterables = [("subject_id", subjects)]
    else:
        infosource.iterables = [("subject_id", [subjects[subjects.index(subject)]]),
                                ("bold_dirs", bolddirs_by_subject[subject])]
    
 
    #setup the datasources
    datasource = nipype.pipeline.engine.Node(
                    nipype.interfaces.io.DataGrabber(infields=["subject_id", "bold_dirs"],
                                                     outfields=["anat", "bold"]),
                                             name="datasource")
    
    datasource.inputs.base_directory = data_directory
    datasource.inputs.template = "*"
    datasource.inputs.field_template = {
                                        "bold" : "%s/bold/%s/bold.nii*"
                                        }
    datasource.inputs.template_args = {
                                       "bold" : [["subject_id", "bold_dirs"]]
                                       }
    datasource.inputs.sort_filelist = True
    
    workflow = nipype.pipeline.engine.Workflow(name="openfmri")
    
    if work_directory is None:
        work_directory = os.path.join(os.getcwd(),'working')
        
    workflow.base_dir = work_directory
    
    workflow.connect(infosource, "subject_id", datasource, "subject_id")
    workflow.connect(infosource, "bold_dirs", datasource, "bold_dirs")

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

def analyze_openfmri_dataset(data_directory, subject=None, model_id=None, task_id=None,  work_directory=None, xnat_datasource=False,
                             run_motion_correction=False, run_skull_strip=False, run_melodic=False, run_autorecon=False):
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
    datasource.inputs.template = "*"
    datasource.inputs.field_template = {'dicoms' : '%s/*dcm'}
    datasource.inputs.template_args = {'dicoms' : [['scan_id']]}
    datasource.inputs.sorted = False
    datasource.inputs.sort_filelist = False

    
    #build the nifti converter
    converter = nipype.pipeline.engine.Node(
                    nipype.interfaces.dcm2nii.Dcm2nii(terminal_output="file"),
                    name="nifticonvert"
    )
    
    #this is a function interface that grabs the dicom info for later
    dicomheaders = nipype.pipeline.engine.Node(
                        nipype.interfaces.utility.Function(input_names=["dicom_file"],
                                                           output_names=["dicom_header"],
                                                           function=get_dicom_headers),
                                               name="dicomheaders")


    openfmri_organization = nipype.pipeline.engine.Node(
                                nipype.interfaces.utility.Function(input_names=["dicom_header","niftis","base_directory"],
                                                                   output_names = ["destination"],
                                                                   function=direct_nifti_to_directory),
                                                        name="openfmri_organization")
    openfmri_organization.inputs.base_directory = os.path.join(openfmri_subject_directory, subject_id)

    #just a test sink for now
    datasink = nipype.pipeline.engine.Node(
                    nipype.interfaces.io.DataSink(base_directory=os.path.join(openfmri_subject_directory, subject_id)),
                    name="datasink")

    datasink.inputs.substitutions = [('_scan_id_', '')] 
    
    workflow = nipype.pipeline.engine.Workflow(name='dicom2nifti')
    workflow.connect(infosource, 'scan_id', datasource, 'scan_id')
    #lambda function that just runs pop() on any list you hand it
    pop_last = lambda x: x.pop() if type(x) is list else x
    workflow.connect(datasource, ('dicoms', pop_last), converter, "source_names")
    workflow.connect(datasource, ('dicoms', pop_last), dicomheaders, "dicom_file")
    workflow.connect(dicomheaders, "dicom_header", openfmri_organization, "dicom_header")
    workflow.connect(converter, "converted_files", openfmri_organization, "niftis")
    workflow.connect(converter, "converted_files", datasink, "niftis")

    workflow.base_dir = os.environ.get("SCRATCH", "/tmp")
    workflow.run()
    

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


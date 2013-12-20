import nipype.interfaces.base
import enthought.traits.api
import os

class DTIQAInputSpec(nipype.interfaces.base.CommandLineInputSpec):
    in_file = nipype.interfaces.base.File(desc = "DTI File", exists = True, mandantory = True, argstr="%s 0.2 .")

class DTIQAOutputSpec(nipype.interfaces.base.TraitedSpec):
    snr = nipype.interfaces.base.File(desc = "Signal PNG", exists = True)
    slice_corr = nipype.interfaces.base.File(desc = "Slice correlation PNG", exists = True)
    interleave_corr = nipype.interfaces.base.File(desc = "Interleave correlation PNG", exists = True)
    fa = nipype.interfaces.base.File(desc = "FA PNG", exists = True)
    fd = nipype.interfaces.base.File(desc = "FD PNG", exists = True)
    sse = nipype.interfaces.base.File(desc = "SSE PNG", exists = True)
    worst_grad = nipype.interfaces.base.File(desc = "Worst Gradient PNG", exists = True)
    report = nipype.interfaces.base.File(desc = "QA Report PDF", exists = True)

class DTIQATask(nipype.interfaces.base.CommandLine):
    input_spec = DTIQAInputSpec
    output_spec = DTIQAOutputSpec
    cmd = "dtiqa.py"

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['snr'] = os.path.abspath('snr.png')
        outputs['slice_corr'] = os.path.abspath('slicecorr.png')
        outputs['interleave_corr'] = os.path.abspath('interleavecorr.png')
        outputs['fa'] = os.path.abspath('FA.png')
        #I have no idea why one is caps and one isn't, but it is
        outputs['fd'] = os.path.abspath('fd.png')
        outputs['sse'] = os.path.abspath('sse.png')
        outputs['worst_grad'] = os.path.abspath('worst_gradient.png')
        outputs['report'] = os.path.abspath('QA_report.pdf')
        return outputs
    

class FMRIQAInputSpec(nipype.interfaces.base.CommandLineInputSpec):
    in_file = nipype.interfaces.base.File(desc="Motion corrected BOLD file", exists = True,
                                          position = 0, mandantory = True, argstr="%s")
    temporal_resolution = enthought.traits.api.Float(desc="Temporal Resolution of data (TR)", position = 1, 
                                                     mandantory = True, argstr="%s")

class FMRIQAOutputSpec(nipype.interfaces.base.TraitedSpec):
    confound = nipype.interfaces.base.File(desc="confound list", exists=True)
    dvars_png = nipype.interfaces.base.File(desc="DVARS PNG", exists=True)
    dvars_list = nipype.interfaces.base.File(desc="DVARS text list", exists=True)
    fd_list = nipype.interfaces.base.File(desc = "FA PNG", exists = True)
    fd_png = nipype.interfaces.base.File(desc = "FD PNG", exists = True)
    mad = nipype.interfaces.base.File(desc = "mad PNG", exists = True)
    maskmean = nipype.interfaces.base.File(desc = "mask mean PNG", exists=True)
    meanpsd = nipype.interfaces.base.File(desc ="mean psd PNG", exists=True)
    qadata = nipype.interfaces.base.File(desc ="qa output CSV", exists=True)
    qa_report = nipype.interfaces.base.File(desc ="Full QA report PDF", exists=True)
    spike = nipype.interfaces.base.File(desc ="spike PNG", exists=True)
    voxcv = nipype.interfaces.base.File(desc ="voxcv PNG", exists=True)
    voxmean = nipype.interfaces.base.File(desc ="vox mean PNG", exists=True)
    voxsfnr_volume = nipype.interfaces.base.File(desc ="vox fnr Nifti", exists=True)
    voxsfnr_png = nipype.interfaces.base.File(desc ="vox fnr PNG", exists=True)
    
class FMRIQATask(nipype.interfaces.base.CommandLine):
    input_spec = FMRIQAInputSpec
    output_spec = FMRIQAOutputSpec
    cmd = "fmriqa.py"
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['confound'] = os.path.abspath('confound.txt')
        outputs['dvars_png'] = os.path.abspath('DVARS.png')
        outputs['dvars_list'] = os.path.abspath('dvars.txt')
        outputs['fd_list'] = os.path.abspath('fd.txt')
        outputs['fd_png'] = os.path.abspath('fd.txt')
        outputs['mad'] = os.path.abspath('mad.png')
        outputs['maskmean'] = os.path.abspath('maskmean.png')
        outputs['meanpsd'] = os.path.abspath('meanpsd.png')
        outputs['qadata'] = os.path.abspath('qadata.csv')
        outputs['spike'] = os.path.abspath('spike.png')
        outputs['voxcv'] = os.path.abspath('voxcv.png')
        outputs['voxmean'] = os.path.abspath('voxmean.png')
        outputs['voxsfnr'] = os.path.abspath('voxsfnr.nii.gz')
        outputs['voxsfnr'] = os.path.abspath('voxsfnr.png')
        return outputs
 
if __name__ == '__main__':
    test = DTIQATask(in_file="DTI_1.nii.gz")
    print test.cmdline
    test.run()
    test = FMRIQATask(in_file="bold_mcf.nii.gz")
    print test.cmdline
    test.run()


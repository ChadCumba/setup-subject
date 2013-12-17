import nipype.interfaces.base
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
 
if __name__ == '__main__':
    test = DTIQATask(in_file="DTI_1.nii.gz")
    print test.cmdline
    test.run()

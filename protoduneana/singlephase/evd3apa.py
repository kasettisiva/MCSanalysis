#!/usr/bin/env python

import h5py
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import sys
import math

def main(argv):

        inputfile = 'r5387_e89634.h5'
        t0 = 0
        t1 = 6000
        w0 = 0
        w1 = 480*3
        sc = 50
        do_save_image = 0

        args = argv[1:]
        while len(args) > 0:
                if args[0] == '-f' and len(args) > 1:
                        inputfile = args[1]
                        del args[0:2]
                elif args[0] == '-w0' and len(args) > 1:
                        w0 = int(args[1])
                        del args[0:2]
                elif args[0] == '-w1' and len(args) > 1:
                        w1 = int(args[1])
                        del args[0:2]
                elif args[0] == '-t0' and len(args) > 1:
                        t0 = int(args[1])
                        del args[0:2]
                elif args[0] == '-t1' and len(args) > 1:
                        t1 = int(args[1])
                        del args[0:2]
                elif args[0] == '-sc' and len(args) > 1:
                        sc = int(args[1])
                        del args[0:2]
                elif args[0] == '-save':
                        do_save_image = 1
                        del args[0]
                else:
                        print('Unknown option %s' % args[0])
                        return 1

        f = h5py.File(inputfile, 'r')

        #fsize = math.log(sc)*5
        #mpl.rcParams.update({'font.size' : fsize})
        mpl.rcParams.update({'font.size' : (w1-w0)/sc/1.6})
        
        #f = h5py.File("r5779_e12360.h5", 'r')
        #f = h5py.File("r5770_e59001.h5", 'r')
        #f = h5py.File("r5826_e83959.h5", 'r')
        #f = h5py.File("r5759_e58084.h5", 'r')
        #f = h5py.File("r5244_e10488.h5", 'r')
        #f = h5py.File("r5809_e10747.h5", 'r')
        
        #f = h5py.File("r5439_e13.h5", 'r')
        #f = h5py.File("r5145_e81569.h5", 'r')
        #f = h5py.File("r5458_e21389.h5", 'r')
        #f = h5py.File("r5432_e50486.h5", 'r')
        #f = h5py.File("r5772_e15132.h5", 'r')
        #f = h5py.File("r5815_e962.h5", 'r')

        #f = h5py.File("r5770_e50648.h5", 'r')
        #f = h5py.File("r5203_e1290.h5", 'r')
        #f = h5py.File("r5759_e84253.h5", 'r')
        #f = h5py.File("r5759_e84081.h5", 'r')
        #f = h5py.File("r5759_e58.h5", 'r')
        #f = h5py.File("r5759_e77.h5", 'r')
        #f = h5py.File("r5145_e27422.h5", 'r')
        #f = h5py.File("r5145_e27543.h5", 'r')
        #f = h5py.File("r5387_e89634.h5", 'r')


        n1 = np.array(f["wiresigs/adc"][:])
        evttime = np.array(f["evtids/evttime"][:])
        run = np.array(f["evtids/eid"][:])
        from datetime import datetime
        timestamp = int(evttime[0])
        dt_object = datetime.fromtimestamp(timestamp)
        import matplotlib
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        y = np.reshape(n1, (n1.shape[0]//6000, 6000))
        z = np.flip(y.transpose(),0)
        #t0 = 0
        #t1 = 5999
        #t0 = 2750
        #t1 = 4999
        #t0 = 2000
        #t1 = 4999
        #w0 = 0
        #w1 = 480*3-1
        #w0 = 0
        #w1 = 250
        z1 = z[6000-t1:6000-t0, w0:w1]
        #mpl.rcParams["image.origin"] = 'lower'
        #im = ax.imshow(z1,  extent=[w0, w1, t0, t1], aspect='auto', vmax=10.0, vmin=-2.5)
        im = ax.imshow(z1,  extent=[w0, w1, t0, t1], aspect='auto', vmax=10.0, vmin=-2.5)
        im.set_cmap('jet')
        divider = make_axes_locatable(plt.gca())
        cax = divider.append_axes("right", "3%", pad="2%")
        cb = plt.colorbar(im, cax=cax)
        
        #cb = fig.colorbar(im, ax=ax, extend='max', extendrect=True, use_gridspec=True)
        cb.set_label("Charge/tick/channel (ke)")
        fig.set_size_inches((w1-w0+1)*0.479/sc, (t1-t0+1)*0.5*0.16/sc)
        #fig.set_size_inches((w1-w0+1)*0.479/30, (t1-t0+1)*0.5*0.16/30)
        #fig.set_size_inches((w1-w0+1)*0.479/20, (t1-t0+1)*0.5*0.16/20)
        ax.set_xlabel("Wire number")
        ax.set_ylabel("Ticks (0.5 $\mu$s)")
        ax.set_title("ProtoDUNE-SP Run {} Event {} @{} UTC".format(run[0][0], run[0][2], dt_object))
        #plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
        plt.tight_layout()
        if do_save_image:
                plt.savefig("R{}_E{}_T1T5T9_w{}_{}_t{}_{}_sc{}.png".format(run[0][0],run[0][2],w0,w1,t0,t1,sc))
                plt.savefig("R{}_E{}_T1T5T9_w{}_{}_t{}_{}_sc{}.pdf".format(run[0][0],run[0][2],w0,w1,t0,t1,sc))
        plt.show()

if __name__ == '__main__':
        sys.exit(main(sys.argv))

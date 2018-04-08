import numpy as np
import nmrglue as ng
import os
import matplotlib.pyplot as plt

subdirs = [x for x in os.listdir(".") if os.path.isdir(x)]

for subdir in subdirs[:]:
    try:
        dic,data = ng.bruker.read(subdir)
        print("Reading subdir: %s" %(subdir))
        udic = ng.bruker.guess_udic(dic,data)
        #~ print udic
        header_string="ndim:%i\n" %(udic["ndim"])
        for i in range(udic["ndim"]):
            for j,k in udic[i].iteritems():
                #~ print "%s:%s" %(j,k)
                header_string=header_string+"%s:%s\n" %(j,k)
        adv_header_string=""
        for j,k in dic.iteritems():
                if type(dic[j]) == type({}):
                        adv_header_string=adv_header_string+"%s:\n" %(j)
                        for l,m in dic[j].iteritems():
                                if type(dic[j][l]) == 'dict':
                                        adv_header_string=adv_header_string+"%s:\n" %(l)
                                        for p,q in dic[j][l].iteritems():
                                                #~ print "%s:%s" %(p,q)
                                                adv_header_string=adv_header_string+"%s:%s\n" %(p,q)
                                else:
                                        #~ print "%s:%s" %(l,m)
                                        adv_header_string=adv_header_string+"%s:%s\n" %(l,m)
                else:
                        adv_header_string=adv_header_string+"%s:%s\n" %(j,k)
        
        np.shape(data)
        if udic["ndim"] == 1:
                print("Saving 1D Data\n")
                np.savetxt("./%s/fid.ts" %(subdir),np.transpose(np.array([np.arange(0,len(data),1)*1./(udic[0]['sw']),data.real,np.zeros(len(data)),data.imag,np.zeros(len(data))])), header = header_string)
        elif udic["ndim"] == 2:
                print("Saving 2D Data for T1 measurement\n")
                #~ plt.plot(np.arange(0,len(data[0]),1)*1./(2*udic[0]['sw']), data[-1].real)
                #~ plt.plot(np.arange(0,len(data[0]),1)*1./(2*udic[0]['sw']), data[-1].imag)
                phase=180/180.*np.pi
                plt.plot(np.cos(phase)*data[-1].real+np.sin(phase)*data[-1].imag)
                plt.plot(-1*np.sin(phase)*data[-1].real+np.cos(phase)*data[-1].imag)
                plt.show()
                plt.clf()
                plt.plot(np.cos(phase)*np.mean(data[:,72:76].real,axis=1)+np.sin(phase)*np.mean(data[:,72:76].imag,axis=1),"o")
                plt.show()
                np.savetxt("./%s/T1_data.txt" %(subdir),np.cos(phase)*np.mean(data[:,72:76].real,axis=1)+np.sin(phase)*np.mean(data[:,72:76].imag,axis=1), header = header_string)
                #~ print len(data)
                #~ np.savetxt("./%s/fid2d.ts" %(subdir),np.transpose(np.array([np.arange(0,len(data),1)*1./(2*udic[0]['sw']),data.real,data.imag])), header = header_string)
        header_file=open("./%s/adv_header.txt" %(subdir), "w")
        header_file.write("%s" %(adv_header_string))
        header_file.close()
    except Exception as e:
        print("Error in subdir %s: %s" %(subdir,e))


#~ fid_data=np.loadtxt("./FID.txt",delimiter=",")
#~ print(fid_data)
#~ 
#~ zeros=np.zeros(len(fid_data))
#~ 
#~ fid_data_dofft=(np.transpose(np.array([(fid_data[:,0]-1)*1e-5,fid_data[:,1],zeros,fid_data[:,2],zeros])))
#~ 
#~ np.savetxt("FID.ts",fid_data_dofft)


#~ np.savetxt("fid.ts",np.transpose(np.array([np.arange(0,len(data),1)*1./(2*dic['acqus']['SW_h']),data.real,np.zeros(len(data)),data.imag,np.zeros(len(data))])))


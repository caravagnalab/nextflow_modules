//
// FORMATTING SUB-WORKFLOW
//

include { CNA_PROCESSING } from '../../modules/cna2CNAqc/main'
include { VCF_PROCESSING } from '../../modules/vcf2CNAqc/main'


workflow FORMATTER {
    take:
        input
    
    main:
        //tmp = input.map{dataset, patient, sample, f, sm -> f}.view()
        tmp = input.map{dataset, patient, sample, f, sm}.
                branch{ f ->
                        VCF: f.getExtension('.vcf')
                        CNA: f.getExtension('.txt')
                }
        tmp.view()
                //.filter(it.name =~ /(\.vcf|1\.fastq)$/)
                //.view()
        // .filter(it.name =~ /(\.vcf|1\.fastq)$/)
        // .view()
        // println tmp.getClass()
        // //println tmp.f.getClass()

        // string = tmp.toString()
        // println string.getClass()
        // println string
        // myFileChannel = Channel.fromPath(tmp)
        // myFileChannel.getExtension()
        //tmp.getExtension().view()
        //tmp.toString().contains('vcf')

        // c = "no"
        // c.println()
        //out = 'no vcf'
        // if (tmp.toString().split('/').contains('vcf')){
        //     out = VCF_PROCESSING(input)
        //  } else{
        //     out = CNA_PROCESSING(input)
        //  }
       //f = tmp.split('.').contains('.vcf')

        // if (file.getExtension() == '.vcf'){
        //    result =  VCF_PROCESSING(file)
        // } else if (current_input.getExtension() == '.rds'){
        //     VCF_PROCESSING() 
        // } else{
        //     CNA_PROCESSING()
        // }

    emit:
        tmp.view()
        //out
        //out
        //out.println()

}
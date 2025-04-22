function [ struct_irasa ] = func_structure_spectra(cell_in, chlab, cond, ctype, f_list)

N = length(cell_in);
nch = length(chlab);
s_list = {'mixd','frac','osci'};

struct_tmp1 = struct('freq',[],'mixd',[],'frac',[],'osci',[],'Freq_lo',[],'Plaw_lo',[],'Freq_hi',[],'Plaw_hi',[]);
struct_tmp = struct();
for ftype = 1:length(f_list)
    struct_tmp.(f_list{ftype}) = struct_tmp1;
    struct_tmp.(f_list{ftype}).cond = cond;
    struct_tmp.(f_list{ftype}).ctype = ctype;
end


struct_irasa = struct();
for ch = 1:nch
    struct_irasa.(chlab{ch}) = struct_tmp;
end

for ch = 1:nch
    for f = 1:length(f_list)
        for s = 1:length(s_list)
            mat_tmp = [];
            for n = 1:N
                tmp = cell_in{n};
                if isempty(mat_tmp)
                    mat_tmp = zeros(length(tmp.(f_list{f}).freq),N);
                end

                mat_tmp(:,n) = tmp.(f_list{f}).(s_list{s})(:,ch);
            end
            struct_irasa.(chlab{ch}).(f_list{f}).(s_list{s}) = mat_tmp;
        end

        mat_lo_tmp = [];
        for n = 1:N
            tmp = cell_in{n};
            if isempty(mat_lo_tmp)
                mat_lo_tmp = zeros(length(tmp.(f_list{f}).Freq_lo),N);
            end

            mat_lo_tmp(:,n) = tmp.(f_list{f}).Plaw_lo(:,ch);
        end
        struct_irasa.(chlab{ch}).(f_list{f}).Plaw_lo = mat_lo_tmp;

        mat_hi_tmp = [];
        for n = 1:N
            tmp = cell_in{n};
            if isempty(mat_hi_tmp)
                mat_hi_tmp = zeros(length(tmp.(f_list{f}).Freq_hi),N);
            end

            mat_hi_tmp(:,n) = tmp.(f_list{f}).Plaw_hi(:,ch);
        end
        struct_irasa.(chlab{ch}).(f_list{f}).Plaw_hi = mat_hi_tmp;

        struct_irasa.(chlab{ch}).(f_list{f}).freq = tmp.(f_list{f}).freq;
        struct_irasa.(chlab{ch}).(f_list{f}).Freq_lo = tmp.(f_list{f}).Freq_lo;
        struct_irasa.(chlab{ch}).(f_list{f}).Freq_hi = tmp.(f_list{f}).Freq_hi;
    end
end


end
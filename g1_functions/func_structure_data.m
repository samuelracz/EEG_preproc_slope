function [ struct_irasa] = func_structure_data(cell_in,chlab,cond,ctype, f_list)

N = length(cell_in);
nch = length(chlab);
table_H = table(cell(N,1),cell(N,1),zeros(N,1), zeros(N,1), zeros(N,1), 'VariableNames', cat(2,{'cond','ctype'},f_list));
struct_irasa = struct();
for ch = 1:nch
    struct_irasa.(['bb_' chlab{ch}]) = table_H;
    struct_irasa.(['lo_' chlab{ch}]) = table_H;
    struct_irasa.(['hi_' chlab{ch}]) = table_H;
    struct_irasa.(['ar_' chlab{ch}]) = table_H;
    struct_irasa.(['ao_' chlab{ch}]) = table_H;
    struct_irasa.(['lar_' chlab{ch}]) = table_H;
    struct_irasa.(['lao_' chlab{ch}]) = table_H;
end

for n = 1:N
    tmp = cell_in{n};
    for f = 1:length(f_list)
        for ch = 1:nch
            % IRASA
            struct_irasa.(['bb_' chlab{ch}]).cond{n} = cond;
            struct_irasa.(['bb_' chlab{ch}]).ctype{n} = ctype;
            struct_irasa.(['bb_' chlab{ch}]).(f_list{f})(n) = tmp.(f_list{f}).Beta_bb(ch);
            struct_irasa.(['lo_' chlab{ch}]).cond{n} = cond;
            struct_irasa.(['lo_' chlab{ch}]).ctype{n} = ctype;
            struct_irasa.(['lo_' chlab{ch}]).(f_list{f})(n) = tmp.(f_list{f}).Beta_lo(ch);
            struct_irasa.(['hi_' chlab{ch}]).cond{n} = cond;
            struct_irasa.(['hi_' chlab{ch}]).ctype{n} = ctype;
            struct_irasa.(['hi_' chlab{ch}]).(f_list{f})(n) = tmp.(f_list{f}).Beta_hi(ch);
            struct_irasa.(['ar_' chlab{ch}]).cond{n} = cond;
            struct_irasa.(['ar_' chlab{ch}]).ctype{n} = ctype;
            struct_irasa.(['ar_' chlab{ch}]).(f_list{f})(n) = tmp.(f_list{f}).Alpha_raw(ch);
            struct_irasa.(['ao_' chlab{ch}]).cond{n} = cond;
            struct_irasa.(['ao_' chlab{ch}]).ctype{n} = ctype;
            struct_irasa.(['ao_' chlab{ch}]).(f_list{f})(n) = tmp.(f_list{f}).Alpha_osci(ch);
            struct_irasa.(['lar_' chlab{ch}]).cond{n} = cond;
            struct_irasa.(['lar_' chlab{ch}]).ctype{n} = ctype;
            struct_irasa.(['lar_' chlab{ch}]).(f_list{f})(n) = tmp.(f_list{f}).logAlpha_raw(ch);
            struct_irasa.(['lao_' chlab{ch}]).cond{n} = cond;
            struct_irasa.(['lao_' chlab{ch}]).ctype{n} = ctype;
            struct_irasa.(['lao_' chlab{ch}]).(f_list{f})(n) = tmp.(f_list{f}).logAlpha_osci(ch);
        end
    end
end

end
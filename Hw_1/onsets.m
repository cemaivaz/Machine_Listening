
function pos = onsets(peaks,tgap,win_length,shift)

pos = 1+(peaks-1).*tgap + floor(shift*win_length/16);

end
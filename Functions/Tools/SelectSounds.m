function SoundsToUse = SelectSounds(Sounds,param)

switch Sounds
    case 'SpeechMusic'
        SoundsToUse = find(param.snd.SpeechMusicLabels);
    case 'Ferrets'
        SoundsToUse = find(param.snd.FerretLabels);
    case 'NoFerrets'
        SoundsToUse = find(~param.snd.FerretLabels);
    case 'All'
        SoundsToUse = 1:param.snd.Nsound;
    case 'Music'
        SoundsToUse = find(param.snd.MusicLabels);
    case 'Speech'
        SoundsToUse = find(param.snd.SpeechLabels);
    case 'Others'
        SoundsToUse = setdiff(1:param.snd.Nsound,[find(param.snd.FerretLabels) find(param.snd.SpeechMusicLabels)]);
  
    otherwise
        if ~isempty(find(strcmp(Sounds,param.snd.Categories),1))
            
            SoundsToUse =  find(param.snd.CategoryLabels == find(strcmp(Sounds,param.snd.Categories)));
            
        elseif isnumeric(Sounds) && max(Sounds) <= length(param.snd.Categories)
            SoundsToUse=[];
            for cat=1:length(Sounds)
                SoundsToUse = [SoundsToUse find(param.snd.CategoryLabels == Sounds(cat))];
            end
            SoundsToUse=unique(SoundsToUse);
        else
            SoundsToUse=[];
        end
end



if isempty(SoundsToUse)
    error('No matching sounds found...')
end
 
end
function res = clip(input,bot,top)
res = max(min(input,top),bot);
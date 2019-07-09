t = read.table("/home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/data/Output/000_Trx/000_Trx_mu_2negative.pts")
t2 = read.table("/home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/data/Output/000_Trx/000_Trx_mu_2positive.pts")


hist(t$V1)
hist(t2$V1)


s = read.table("/home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/data/Output/000_Trx/000_Trx_mu_100negative.pts")
s2 = read.table("/home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/data/Output/000_Trx/000_Trx_mu_100positive.pts")


hist(s$V1)
hist(s2$V1)

